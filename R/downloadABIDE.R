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
      FALSE
    } else {
      TRUE
    }
  } else {
    TRUE
  }
}

#' Download ABIDE data for analysis.
#'
#' @param outdir A directory to save the output files.
#' @param derivatives Objects to download from the data respository.
#' If freesurfer is selected it downloads full output directories.
#' All derivatives: alff | degree_binarize | degree_weighted | dual_regression |
#'   eigenvector_binarize | eigenvector_weighted | falff | func_mask |
#'   func_mean | func_preproc | lfcd | reho | rois_aal | rois_cc200 |
#'   rois_cc400 | rois_dosenbach160 | rois_ez | rois_ho | rois_tt | vmhc | freesurfer
#' @param pipelines What processing pipeline to use.  All pipelines: ccs | cpac | dparsf | niak
#' @param strategies What processing strategy to use. All strategies: filt_global | filt_noglobal | nofilt_global |
#'   nofilt_noglobal
#' @param force Always re-download existing files.
#' @param QCremove Remove using quality checks.
#' @param ... Arguments passed to mclapply
#' @return Returns a data frame of ABIDE demographic data and paths to the files.
#' @importFrom utils read.csv
#' @importFrom parallel mcmapply
#' @importFrom stats reshape
#' @export
ABIDE <- function(outdir, derivatives=c("alff", "func_mask"),
pipelines='cpac', strategies="filt_global", force = FALSE, QCremove=FALSE, ...) {
  derivatives = unique(derivatives)
  # To manage URL locations for different data types.
  funcDerivatives = c("alff", "degree_binarize", "degree_weighted", "dual_regression", "eigenvector_binarize", "eigenvector_weighted", "falff", "func_mask", "
                func_mean", "func_preproc", "lfcd", "reho", "rois_aal", "rois_cc200", "rois_cc400", "rois_dosenbach160", "rois_ez", "rois_ho", "rois_tt", "vmhc")
  labeldir = file.path('label', c('BA.ctab', 'aparc.annot.a2009s.ctab', 'aparc.annot.ctab', 'lh.BA.annot', 'lh.BA1.label', 'lh.BA2.label', 'lh.BA3a.label', 'lh.BA3b.label', 'lh.BA44.label', 'lh.BA45.label', 'lh.BA4a.label', 'lh.BA4p.label', 'lh.BA6.label', 'lh.MT.label', 'lh.V1.label', 'lh.V2.label', 'lh.aparc.a2009s.annot', 'lh.aparc.annot', 'lh.cortex.label', 'lh.entorhinal_exvivo.label', 'rh.BA.annot', 'rh.BA1.label', 'rh.BA2.label', 'rh.BA3a.label', 'rh.BA3b.label', 'rh.BA44.label', 'rh.BA45.label', 'rh.BA4a.label', 'rh.BA4p.label', 'rh.BA6.label', 'rh.MT.label', 'rh.V1.label', 'rh.V2.label', 'rh.aparc.a2009s.annot', 'rh.aparc.annot', 'rh.cortex.label', 'rh.entorhinal_exvivo.label'))
  mridir = file.path('mri', c('.xdebug_mris_calc', 'T1.mgz', 'aparc+aseg.mgz', 'aparc.a2009s+aseg.mgz', 'aseg.auto.mgz', 'aseg.auto_noCCseg.label_intensities.txt', 'aseg.auto_noCCseg.mgz', 'aseg.mgz', 'brain.finalsurfs.mgz', 'brain.mgz', 'brain_loose.mgz', 'brain_tight.mgz', 'brainmask.auto.mgz', 'brainmask.edit.mgz', 'brainmask.fsinit.mgz', 'brainmask.gcuts.mgz', 'brainmask.loose.mgz', 'brainmask.mgz', 'brainmask.tight.mgz', 'ctrl_pts.mgz', 'filled.mgz', 'lh.ribbon.mgz', 'mri_nu_correct.mni.log', 'norm.mgz', 'nu.mgz', 'nu_noneck.mgz', 'orig.mgz', 'rawavg.mgz', 'rh.ribbon.mgz', 'ribbon.mgz', 'segment.dat', 'talairach.label_intensities.txt', 'talairach.log', 'talairach_with_skull.log', 'transforms/cc_up.lta', 'transforms/talairach.auto.xfm', 'transforms/talairach.lta', 'transforms/talairach.m3z', 'transforms/talairach.m3z.inv.x.mgz', 'transforms/talairach.m3z.inv.y.mgz', 'transforms/talairach.m3z.inv.z.mgz', 'transforms/talairach.xfm', 'transforms/talairach_avi.log', 'transforms/talairach_with_skull.lta', 'transforms/talsrcimg_to_711-2C_as_mni_average_305_t4_vox2vox.txt', 'wm.asegedit.mgz', 'wm.mgz', 'wm.seg.mgz', 'wmparc.mgz'))
  scriptsdir = file.path('scripts',c('TR.list', 'build-stamp.txt', 'ccs_01_anatpreproc.log', 'ccs_01_anatsurfrecon.log', 'ccs_01_funcpreproc.log', 'ccs_02_anatregister.log', 'ccs_02_funcregister.log', 'ccs_03_funcsegment.log', 'ccs_04_funcnuisance.log', 'ccs_05_funcpreproc_final.log', 'ccs_anatproc.sh', 'ccs_funcproc.sh', 'csurfdir', 'lfcd_x.sh', 'ponscc.cut.log', 'recon-all-status.log', 'recon-all.cmd', 'recon-all.done', 'recon-all.env', 'recon-all.env.bak', 'recon-all.local-copy', 'recon-all.log'))
  statsdir = file.path('stats', c('aseg.stats', 'lh.BA.stats', 'lh.aparc.a2009s.stats', 'lh.aparc.stats', 'lh.entorhinal_exvivo.stats', 'rh.BA.stats', 'rh.aparc.a2009s.stats', 'rh.aparc.stats', 'rh.entorhinal_exvivo.stats', 'wmparc.stats'))
  surfdir = file.path('surf', c('lh.area', 'lh.area.fsaverage.mgh', 'lh.area.fsaverage.mris_preproc.log', 'lh.area.fwhm0.fsaverage.mgh', 'lh.area.fwhm10.fsaverage.mgh', 'lh.area.fwhm15.fsaverage.mgh', 'lh.area.fwhm20.fsaverage.mgh', 'lh.area.fwhm25.fsaverage.mgh', 'lh.area.fwhm5.fsaverage.mgh',
                                'lh.area.mid', 'lh.area.pial', 'lh.area.pial.fsaverage.mgh', 'lh.area.pial.fsaverage.mris_preproc.log', 'lh.area.pial.fwhm0.fsaverage.mgh', 'lh.area.pial.fwhm10.fsaverage.mgh', 'lh.area.pial.fwhm15.fsaverage.mgh', 'lh.area.pial.fwhm20.fsaverage.mgh', 'lh.area.pial.fwhm25.fsaverage.mgh',
                                'lh.area.pial.fwhm5.fsaverage.mgh', 'lh.avg_curv', 'lh.curv', 'lh.curv.fsaverage.mgh', 'lh.curv.fsaverage.mris_preproc.log', 'lh.curv.fwhm0.fsaverage.mgh', 'lh.curv.fwhm10.fsaverage.mgh', 'lh.curv.fwhm15.fsaverage.mgh', 'lh.curv.fwhm20.fsaverage.mgh', 'lh.curv.fwhm25.fsaverage.mgh',
                                'lh.curv.fwhm5.fsaverage.mgh', 'lh.curv.pial', 'lh.defect_borders', 'lh.defect_chull', 'lh.defect_labels', 'lh.inflated', 'lh.inflated.H', 'lh.inflated.K', 'lh.inflated.nofix', 'lh.jacobian_white', 'lh.jacobian_white.fsaverage.mgh', 'lh.jacobian_white.fsaverage.mris_preproc.log',
                                'lh.jacobian_white.fwhm0.fsaverage.mgh', 'lh.jacobian_white.fwhm10.fsaverage.mgh', 'lh.jacobian_white.fwhm15.fsaverage.mgh', 'lh.jacobian_white.fwhm20.fsaverage.mgh', 'lh.jacobian_white.fwhm25.fsaverage.mgh', 'lh.jacobian_white.fwhm5.fsaverage.mgh', 'lh.orig', 'lh.orig.nofix',
                                'lh.pial', 'lh.qsphere.nofix', 'lh.smoothwm', 'lh.smoothwm.nofix', 'lh.sphere', 'lh.sphere.reg', 'lh.sulc', 'lh.sulc.fsaverage.mgh', 'lh.sulc.fsaverage.mris_preproc.log', 'lh.sulc.fwhm0.fsaverage.mgh', 'lh.sulc.fwhm10.fsaverage.mgh', 'lh.sulc.fwhm15.fsaverage.mgh',
                                'lh.sulc.fwhm20.fsaverage.mgh', 'lh.sulc.fwhm25.fsaverage.mgh', 'lh.sulc.fwhm5.fsaverage.mgh', 'lh.thickness', 'lh.thickness.fsaverage.mgh', 'lh.thickness.fsaverage.mris_preproc.log', 'lh.thickness.fwhm0.fsaverage.mgh', 'lh.thickness.fwhm10.fsaverage.mgh',
                                'lh.thickness.fwhm15.fsaverage.mgh', 'lh.thickness.fwhm20.fsaverage.mgh', 'lh.thickness.fwhm25.fsaverage.mgh', 'lh.thickness.fwhm5.fsaverage.mgh', 'lh.volume', 'lh.volume.fsaverage.mgh', 'lh.volume.fsaverage.mris_preproc.log', 'lh.volume.fwhm0.fsaverage.mgh',
                                'lh.volume.fwhm10.fsaverage.mgh', 'lh.volume.fwhm15.fsaverage.mgh', 'lh.volume.fwhm20.fsaverage.mgh', 'lh.volume.fwhm25.fsaverage.mgh', 'lh.volume.fwhm5.fsaverage.mgh', 'lh.white', 'rh.area', 'rh.area.fsaverage.mgh', 'rh.area.fsaverage.mris_preproc.log',
                                'rh.area.fwhm0.fsaverage.mgh', 'rh.area.fwhm10.fsaverage.mgh', 'rh.area.fwhm15.fsaverage.mgh', 'rh.area.fwhm20.fsaverage.mgh', 'rh.area.fwhm25.fsaverage.mgh', 'rh.area.fwhm5.fsaverage.mgh', 'rh.area.mid', 'rh.area.pial', 'rh.area.pial.fsaverage.mgh',
                                'rh.area.pial.fsaverage.mris_preproc.log', 'rh.area.pial.fwhm0.fsaverage.mgh', 'rh.area.pial.fwhm10.fsaverage.mgh', 'rh.area.pial.fwhm15.fsaverage.mgh', 'rh.area.pial.fwhm20.fsaverage.mgh', 'rh.area.pial.fwhm25.fsaverage.mgh', 'rh.area.pial.fwhm5.fsaverage.mgh',
                                'rh.avg_curv', 'rh.curv', 'rh.curv.fsaverage.mgh', 'rh.curv.fsaverage.mris_preproc.log', 'rh.curv.fwhm0.fsaverage.mgh', 'rh.curv.fwhm10.fsaverage.mgh', 'rh.curv.fwhm15.fsaverage.mgh', 'rh.curv.fwhm20.fsaverage.mgh', 'rh.curv.fwhm25.fsaverage.mgh', 'rh.curv.fwhm5.fsaverage.mgh',
                                'rh.curv.pial', 'rh.defect_borders', 'rh.defect_chull', 'rh.defect_labels', 'rh.inflated', 'rh.inflated.H', 'rh.inflated.K', 'rh.inflated.nofix', 'rh.jacobian_white', 'rh.jacobian_white.fsaverage.mgh', 'rh.jacobian_white.fsaverage.mris_preproc.log', 'rh.jacobian_white.fwhm0.fsaverage.mgh',
                                'rh.jacobian_white.fwhm10.fsaverage.mgh', 'rh.jacobian_white.fwhm15.fsaverage.mgh', 'rh.jacobian_white.fwhm20.fsaverage.mgh', 'rh.jacobian_white.fwhm25.fsaverage.mgh', 'rh.jacobian_white.fwhm5.fsaverage.mgh', 'rh.orig', 'rh.orig.nofix', 'rh.pial', 'rh.qsphere.nofix', 'rh.smoothwm',
                                'rh.smoothwm.nofix', 'rh.sphere', 'rh.sphere.reg', 'rh.sulc', 'rh.sulc.fsaverage.mgh', 'rh.sulc.fsaverage.mris_preproc.log', 'rh.sulc.fwhm0.fsaverage.mgh', 'rh.sulc.fwhm10.fsaverage.mgh', 'rh.sulc.fwhm15.fsaverage.mgh', 'rh.sulc.fwhm20.fsaverage.mgh', 'rh.sulc.fwhm25.fsaverage.mgh',
                                'rh.sulc.fwhm5.fsaverage.mgh', 'rh.thickness', 'rh.thickness.fsaverage.mgh', 'rh.thickness.fsaverage.mris_preproc.log', 'rh.thickness.fwhm0.fsaverage.mgh', 'rh.thickness.fwhm10.fsaverage.mgh', 'rh.thickness.fwhm15.fsaverage.mgh', 'rh.thickness.fwhm20.fsaverage.mgh', 'rh.thickness.fwhm25.fsaverage.mgh',
                                'rh.thickness.fwhm5.fsaverage.mgh', 'rh.volume', 'rh.volume.fsaverage.mgh', 'rh.volume.fsaverage.mris_preproc.log', 'rh.volume.fwhm0.fsaverage.mgh', 'rh.volume.fwhm10.fsaverage.mgh', 'rh.volume.fwhm15.fsaverage.mgh', 'rh.volume.fwhm20.fsaverage.mgh', 'rh.volume.fwhm25.fsaverage.mgh', 'rh.volume.fwhm5.fsaverage.mgh', 'rh.white'))
  fsfiles = c(labeldir, mridir, scriptsdir, statsdir,surfdir)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # download demographic info
  demodir <- file.path(outdir, "demographic")
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
  meta = subset(meta, file_id != "no_filename")
  if(QCremove){
    meta <- subset(meta,
                    (
                     qc_func_rater_2 != "fail" | qc_anat_rater_2 != "fail" |
                       qc_anat_rater_3 != "fail" | qc_rater_1 != "fail"
                   ))
  }

  # Formatting specific to abide
  meta$sex <- factor(meta$sex)
  meta$dx_group <- factor(meta$dx_group, labels = c("asd", "hc"))

  # create directory for imaging data
  neurodir <- file.path(outdir, "neuroimaging")
  if (!dir.exists(neurodir)) {
    dir.create(neurodir)
  }

  # copy template file if needed
  templatefile <- file.path(neurodir, "MNI152_T1_3mm.nii.gz")
  if (!exists(templatefile)) {
    srcfile <- file.path(find.package("PCP"), "extdata", "MNI152_T1_3mm.nii.gz")
    file.copy(srcfile, templatefile, overwrite = TRUE)
  }

  # compile list of files to fetch
  baseurl <- "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs"
  # For functional data
  if(all(derivatives %in% funcDerivatives)){
  files <- expand.grid(pipeline = pipelines, strategy = strategies,
                       derivative = derivatives, file_id = meta$file_id)
  files$destdir <- file.path(neurodir, files$pipeline, files$strategy, files$derivative)
  # Use .nii.gz for all derivatives except for the ROI time
  # series files, which end in .1D (these derivative names begin with rois_).
  files$filename <- paste0(files$file_id, "_", files$derivative, ifelse(grepl('roi', files$derivative), ".1D", '.nii.gz'))
  files$destfile <- file.path(files$destdir, files$filename)
  files$url <- file.path(baseurl, files$pipeline, files$strategy, files$derivative, files$filename)
  # for freesurfer data
  } else if(all(derivatives %in% 'freesurfer')){
    files <- expand.grid(pipeline = 'freesurfer', strategy = '5.1',
                         derivative = fsfiles, file_id = meta$file_id)
    files$filename <- files$derivative
    files$destfile <- file.path(neurodir, files$pipeline, files$strategy, files$file_id, files$filename)
    files$destdir <- dirname(files$destfile)
    files$url <- file.path(baseurl, files$pipeline, files$strategy, files$file_id, files$filename)

  } else {
    stop('Functional and structural derivatives must be downloaded separately.')
  }

  # create destination directories if needed
  destdirs <- unique(files$destdir)
  for (i in which(!dir.exists(destdirs))) {
    dir.create(destdirs[i], recursive = TRUE)
  }

  # download files if needed
  files$success = unlist(mcmapply(downloadABIDEFile, url=files$url, destfile=files$destfile, MoreArgs=c('force'=force, list(...) ) ))

  rfiles = reshape(files, v.names = c('destdir', 'url', 'destfile', 'filename', 'success'),
                   timevar='derivative', idvar='file_id', direction='wide', sep='_')

  # merge with meta
  meta = merge(meta, rfiles)
  return(meta)
}
