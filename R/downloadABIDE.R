#' Download ABIDE data for analysis.
#'
#' @param url Undocumented arguments from Jeremy.
#' @param destfile Undocumented arguments from Jeremy.
#' @param force Undocumented arguments from Jeremy.
#' @return Returns a file?
#' @importFrom httr write_disk
#'
#' @examples
#'
#' # downloadABIDE('./', derivatives='falff')
downloadABIDEFile <- function(url, destfile, force) {
  fetch <- force || !file.exists(destfile)
  if (fetch) {
    message("Downloading ", destfile, "...")
    res <- httr::GET(url, write_disk(destfile, overwrite = TRUE))
    if (res$status != 200) {
      stop("Failed to fetch url: ", url)
    }
  }
}

#' Download ABIDE data for analysis.
#'
#' @param outdir A directory to save the output files.
#' @param derivatives Objects to download from the data respository. All derivatives: alff | degree_binarize | degree_weighted | dual_regression |
#'   eigenvector_binarize | eigenvector_weighted | falff | func_mask |
#'   func_mean | func_preproc | lfcd | reho | rois_aal | rois_cc200 |
#'   rois_cc400 | rois_dosenbach160 | rois_ez | rois_ho | rois_tt | vmhc
#' @param pipelines What processing pipeline to use.  All pipelines: ccs | cpac | dparsf | niak
#' @param strategies What processing strategy to use. All strategies: filt_global | filt_noglobal | nofilt_global |
#'   nofilt_noglobal
#' @param force Always re-download existing files.
#' @return Returns a data frame of ABIDE demographic data.
#' @importFrom utils read.csv
#' @export
downloadABIDE <- function(outdir, derivatives=c("alff", "func_mask"),
pipelines='cpac', strategies="filt_global", force = FALSE) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  oldwd <- setwd(outdir)
  on.exit(setwd(oldwd), add = TRUE)

  # download demographic info
  demodir <- file.path("abide", "demographic")
  if (!dir.exists(demodir)) {
    dir.create(demodir, recursive = TRUE)
  }

  # See http://preprocessed-connectomes-project.org/abide/download.html
  metafile <- file.path(demodir, "Phenotypic_V1_0b_preprocessed1.csv")
  metaurl <- "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Phenotypic_V1_0b_preprocessed1.csv"
  downloadABIDEFile(metaurl, metafile, force)

  meta <- read.csv(metafile, na.strings = c("-9999", ""))
  names(meta) <- tolower(names(meta))

  # exclude records without files and/or with poor quality control ratings
  meta <- subset(meta,
    file_id != "no_filename" & (
      qc_func_rater_2 != "fail" | qc_anat_rater_2 != "fail" |
      qc_anat_rater_3 != "fail" | qc_rater_1 != "fail"
    ))

  meta$sex <- factor(meta$sex)
  meta$dx_group <- factor(meta$dx_group, labels = c("asd", "hc"))

  # create directory for imaging data
  neurodir <- file.path("abide", "neuroimaging")
  if (!dir.exists(neurodir)) {
    dir.create(neurodir)
  }

  # copy template file if needed
  templatefile <- file.path(neurodir, "MNI152_T1_3mm.nii.gz")
  if (!exists(templatefile)) {
    srcfile <- file.path(find.package("NIsim"), "extdata", "MNI152_T1_3mm.nii.gz")
    file.copy(srcfile, templatefile, overwrite = TRUE)
  }

  # compile list of files to fetch
  # NOTE: excluding strategy from destination directory structure since there's
  #       only 1
  baseurl <- "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs"
  files <- expand.grid(pipeline = pipelines, strategy = strategies,
                       derivative = derivatives, file_id = meta$file_id)
  files$destdir <- file.path(neurodir, files$pipeline, files$derivative)
  files$filename <- paste0(files$file_id, "_", files$derivative, ".nii.gz")
  files$destfile <- file.path(files$destdir, files$filename)
  files$url <- paste(baseurl, files$pipeline, files$strategy, files$derivative, files$filename, sep = "/")

  # create destination directories if needed
  destdirs <- unique(files$destdir)
  for (i in which(!dir.exists(destdirs))) {
    dir.create(destdirs[i], recursive = TRUE)
  }

  # download files if needed
  for (i in 1:nrow(files)) {
    downloadABIDEFile(files$url[i], files$destfile[i], force)
  }

  return(meta)
}
