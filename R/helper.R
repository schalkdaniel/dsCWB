#' @title Check if symbol is present in global envir.
#' @param symbol (`character(1L)`)\cr
#'   Symbol specifying something.
#' @export
checkSymbol = function(symbol) {
  checkmate::assertCharacter(symbol, len = 1L, any.missing = FALSE)
  fss = ls(envir = .GlobalEnv)
  e = try(eval(parse(text = symbol), envir = .GlobalEnv), silent = TRUE)
  if (inherits(e, "try-error"))
    stop("Cannot find symbol ", symbol, " in .GlobalEnv")
}


#' @title Paste a numeric vector into a character vector
#' @description
#'   Converts a numerical vector into one string with the values separated
#'   by the character specified in `sep` (default = /).
#' @param x (`numeric()`)\cr
#'   Numeric vector converted into the string.
#' @param sep (`character(1L)`)\cr
#'   Character used to separate numerical values.
#' @param digits (`integer(1L)`)\cr
#'   Significant digits used in the string.
#' @examples
#' x = rnorm(20)
#' xchar = parseParams(x)
#' @export
parseParams = function(x, sep = "/", digits = 22L) {
  checkmate::assertNumeric(x = x, any.missing = FALSE)
  checkmate::assertCharacter(x = sep, len = 1L, any.missing = FALSE)
  checkmate::assertCount(x = digits, positive = TRUE)

  xc = format(x, digits = digits)
  if (! is.null(names(x))) {
    xc = paste0(names(x), "=", xc)
  }
  xchar = gsub(" ", "", paste(xc, collapse = sep), perl = TRUE)
  return(xchar)
}

#' @title Deparse a character vector representation
#' @description
#'   Converts a string created by `parseParams` into the corresponding
#'   numeric vector.
#' @param xchar (`character(1L)`)\cr
#'   Character string obtained by `parseParams`.
#' @param sep (`character(1L)`)\cr
#'   Character used to separate numerical values.
#' @examples
#' x = rnorm(20)
#' xchar = parseParams(x)
#' xnew = deparseParams(xchar)
#' all.equal(x, xnew)
#' @export
deparseParams = function(xchar, sep = "/") {
  checkmate::assertCharacter(x = xchar, len = 1L, any.missing = FALSE)
  checkmate::assertCharacter(x = sep, len = 1L, any.missing = FALSE)

  if (! grepl(sep, xchar)) stop("Cannot find separator \"", sep, "\" in xchar.")

  xsplit = strsplit(xchar, sep)[[1]]
  cl = paste0("c(", paste(xsplit, collapse = ", "), ")")
  return(eval(parse(text = cl)))
}

#' @title Serialize R object
#' @description This function serializes a given R object and creates a character string
#'   containing the binary of the object. This object can be send to the DataSHIELD servers.
#' @param obj (arbitrary R object) Object which should be send to DataSHIELD.
#' @param obj_name (`character(1L)`) Name of the object (default is `NULL`). If name is set to
#'   `NULL`, then the object name passed to the function is used.
#' @param check_serialization (`logical(1L)`) Check if the serialized model can be deserialized
#'   locally (default is `TRUE`).
#' @return Character of length 1 containing the serialized object as string.
#' @author Daniel S.
#' @examples
#' mod = lm(Sepal.Width ~ ., data = iris)
#' bin = encodeObject(mod)
#' substr(bin, 1, 50)
#' @export
encodeObject = function(obj, obj_name = NULL, check_serialization = TRUE) {
  checkmate::assertCharacter(obj_name, len = 1L, null.ok = TRUE, any.missing = FALSE)
  checkmate::assertLogical(check_serialization, len = 1L)

  if (is.null(obj_name)) obj_name = deparse(substitute(obj))

  obj_binary = serialize(obj, connection = NULL)
  obj_binary_str = as.character(obj_binary)
  obj_binary_str_collapsed = paste(obj_binary_str, collapse = "")

  ## Pre check if object serialization works locally:
  if (check_serialization) {
    # get object back from serialization
    obj_b = decodeBinary(obj_binary_str_collapsed)
    if (! all.equal(obj, obj_b)) stop("Model cannot serialized and deserialized into equal object!")
  }

  osize = utils::object.size(obj_binary_str_collapsed) / 1024^2
  if (osize > 1) {
    message("[", Sys.time(), "] Your object is bigger than 1 MB (", round(osize, 1),
      " MB). Uploading larger objects may take some time.")
  }
  names(obj_binary_str_collapsed) = obj_name

  return(obj_binary_str_collapsed)
}

#' @title Deserialize object
#' @description Decode a given string of a serialized object.
#' @param bin (`character(1L)`) Binary string value containing the serialized model.
#' @param package (`character(1L)`) Package required for object deserialization (default is `NULL`).
#' @return Deserialized object from `bin`
#' @author Daniel S.
#' @examples
#' mod = lm(Sepal.Width ~ ., data = iris)
#' bin = encodeObject(mod)
#' mod_b = decodeBinary(bin)
#' all.equal(mod, mod_b)
#' @export
decodeBinary = function(bin, package = NULL) {
  checkmate::assertCharacter(bin, len = 1L, null.ok = FALSE, any.missing = FALSE)
  checkmate::assertCharacter(package, len = 1L, null.ok = TRUE, any.missing = FALSE)

  # Check if model is installed and install if not:
  if (! is.null(package) && require(package, quietly = TRUE, character.only = TRUE)) {
    stop("Package '", package, "' is not installed. Please install it or contact the administrator to do this for you.")
  }

  binary_str_deparse = substring(bin, seq(1, nchar(bin), 2), seq(2, nchar(bin), 2))

  raw = as.raw(as.hexmode(binary_str_deparse))
  obj = unserialize(raw)

  return(obj)
}

#' @title Get the session information of the DataSHIELD server
#' @description This method returns `sessionInfo()` from the used DataSHIELD servers.
#'   The main purpose is for testing and checking the environment used on the remote servers.
#' @return list of session infos returned from `sessionInfo()` of each machine
#' @author Daniel S.
#' @export
getDataSHIELDInfo = function() {
  out = list(
    session_info = utils::sessionInfo(),
    pcks = utils::installed.packages())

  return(out)
}

# Get `datashield.privacyLevel` from DESCRIPTION file. Note that we do not set the option
# as DataSHIELD does because of the risk of highjacking the R environment. Instead, when
# a function is called that uses the privacy level, the function gets it directly from the
# DESCRIPTION file.
.getPrivacyLevel = function() {
  pl = utils::packageDescription("dsBinVal")$Options
  pl = as.integer(gsub("\\D", "", pl))
  if (is.na(pl)) stop("No privacy level specified in DESCRIPTION.")
  return(pl)
}

.tryOPALConnection = function(expr) {
  conns = suppressMessages(try(expr, silent = TRUE))
  if (inherits(conns, "opal")) {
    return(conns)
  } else {
    return("Was not able to establish connection")
  }
}
