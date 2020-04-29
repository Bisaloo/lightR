#' Extract spectral data from spectra files
#'
#' Finds and imports reflectance/transmittance/absorbance data from spectra
#' files in a given location.
#'
#' @inheritParams lr_parse_generic
#'
#' @param where Folder in which files are located (defaults to current working
#' directory).
#' @param ext File extension to be searched for, without the "." (defaults to
#' `txt`). You can also use a character vector to specify multiple file
#' extensions.
#' @param lim A vector with two numbers determining the wavelength limits to be
#' considered (defaults to `c(300, 700)`).
#' @param subdir Should subdirectories within the `where` folder be included in
#' the search? (defaults to `FALSE`).
#' @param subdir.names Should subdirectory path be included in the name of the
#' spectra? (defaults to `FALSE`).
#' @param cores deprecated. See [future::plan()] for more details on how to
#' customise your parallelisation strategy.
#' @param ignore.case Should the extension search be case insensitive? (defaults
#' to `TRUE`)
#' @param interpolate Boolean indicated whether spectral data should be
#' interpolated and pruned at every nanometre. Note that this option can only
#' work if all input data samples the same wavelengths. Defaults to `TRUE`.
#'
#' @details
#' You can customise the type of parallel processing used by this function with
#' the [future::plan()] function. This works on all operating systems, as well
#' as high performance computing (HPC) environment. Similarly, you can customise
#' the way progress is shown with the [progressr::handlers()] functions
#' (progress bar, acoustic feedback, nothing, etc.)
#'
#' @return A data.frame, containing the wavelengths in the first column and
#' individual imported spectral files in the subsequent columns.
#'
#' @export
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom stats approx
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#'
#' @seealso [pavo::getspec()]
#'
#' @examples
#' lr_get_spec(system.file("testdata", package = "lightr"), ext = "jdx")
#'
lr_get_spec <- function(where = getwd(), ext = "txt", lim = c(300, 700),
                        decimal = ".", sep = NULL, subdir = FALSE,
                        subdir.names = FALSE, cores = NULL, ignore.case = TRUE,
                        interpolate = TRUE) {

  if (!missing(cores)) {
    warning("'cores' argument is deprecated. See ?future::plan for more info ",
            "about how you can choose your parallelisation strategy.",
            call. = FALSE)
  }

  extension <- paste0("\\.", ext, "$", collapse = "|")

  file_names <- list.files(where,
    pattern = extension, ignore.case = ignore.case,
    recursive = subdir, include.dirs = subdir
  )

  # This step is needed to ensure reproducibility between locales and platforms
  file_names <- sort(file_names, method = "radix")

  nb_files <- length(file_names)

  if (nb_files == 0) {
    warning('No files found. Try a different value for argument "ext".',
            call. = FALSE)
    return(NULL)
  }

  files <- file.path(where, file_names)

  if (!subdir.names) {
    file_names <- basename(file_names)
  }

  specnames <- file_path_sans_ext(file_names)

  range <- seq(lim[1], lim[2])

  message(nb_files, " files found; importing spectra:")

  if (!interpolate) {
    gsp <- function(f) {
      dispatch_parser(f, decimal = decimal, sep = sep)[[1]]
    }
  } else {
    gsp <- function(f) {
      df <- dispatch_parser(f, decimal = decimal, sep = sep)[[1]]

      # Trim now because otherwise, approx() can fill the region of interest
      # with bogus data (e.g., if the data is complete between 200-300nm and
      # 800-1200nm and your region of interested if 300-700 nm).
      bounds <- which(df$wl >= lim[1] & df$wl <= lim[2])

      if (length(bounds) == 0) {
        warning(f, " does not contain spectral data over the provided wl range",
                call. = FALSE)
        return(NULL)
      }

      df <- df[c(min(bounds)-1, bounds, max(bounds)+1), ]

      approx(df[, "wl"], df[, "processed"],
             xout = range, ties = "ordered")$y
    }
  }

  with_progress({
  p <- progressor(along = files)
  tmp <- future_lapply(files, function(x) {
    p()
    tryCatch(gsp(x),
             error = function(e) NULL)
    })
  })

  whichfailed <- which(vapply(tmp, is.null, logical(1)))

  if (length(whichfailed) == nb_files) {
    warning("File import failed.\n",
            "Check input files and function arguments.", call. = FALSE)
    return(NULL)
  } else if (length(whichfailed) > 0) {

    warning("Could not import one or more files:\n",
      paste0(files[whichfailed], "\n"),
      call. = FALSE
    )

    specnames <- specnames[-whichfailed]
  }

  if (interpolate) {
    # We convert to data.frame before adding the wl to preserve wl class (int or
    # altrep)
    tmp <- list2DF(tmp)

    final <- cbind(range, tmp)
  } else {
    if (length(unique(lapply(tmp, function(x) x[, "wl"])))!=1) {
      stop("'interpolate = FALSE' can only work if all input files sample the ",
           "same wavelengths.", call. = FALSE)
    }
    final <- list2DF(lapply(tmp, function(x) x[, "processed"]))

    final <- cbind(tmp[[1]][, "wl"], final)
  }

  colnames(final) <- c("wl", specnames)

  final <- final[final$wl <= lim[2] & final$wl >= lim[1],]

  # TODO: should uninterpolated spectra be marked as rspec objects?
  class(final) <- c("rspec", "data.frame")

  return(final)
}
