###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
#      name:  clean_copre_folder.R                                             #
# timestamp:  2026-07-30T11:43Z                                                #
#    author:  Tillman Reuter (orcid.org/0009-0002-1156-9983)                   #
#   license:  GPLv3 (see LICENSE.md)                                           #
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%

# READ ME ######################################################################
'
 PURPOSE:
   This script reduces the CoPre2 output folder to just essential files. These
   can be manually specified or instead are the following:
    - trajectory
    - SLAM-report
    - 2 images
    - point cloud .las
  There is a function aimed at the Result folder

 DEPENDENCIES:
   none

 USAGE:
   Source this file in your R session, then call `clean_copre_folder()`.
'

###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
# LOAD AND INSTALL PACKAGES ####################################################
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%

# no dependencies

###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
# FUNCTIONS ####################################################################
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%

#' Retain from a CoPre2 scan folder only the relevant files
#'
#' This function performs a standardized filtering of files and subdirectories
#' of the CoPre2 scan folder, deletes irrelevant ones and collapses the structure.
#'
#' @param folder   Character. Path to the CoPre2 scan folder (something like "C:/.../CHC.../2026-...")
#' @param retain   Character (or Character vector). Regex patterns of relevant
#'                 files to keep. Defaults to
#'                 c("^tracks\\.txt$",
#'                 "^Accuracy_XY\\.jpg$",
#'                 "^[0-9]{8}.*\\.las$",
#'                 "SlamReport\\.html$")
#'
#' @return Returns NULL The actual output is the manipulation of the files and
#'         directories on disk.
#'
#' @export
#' @examples
#' \dontrun{
#' clean_copre_folder("C:/MLS/costa_rica/CHC_20261021111716/2026-11-23-123456(SN-005C)")
#' }
clean_copre_scan <- function(folder, retain=NULL) {

  # Sanity checks
  if (is.null(folder)) stop('Must provide a folder name (something like "C:/.../CHC.../2026-...")')
  if (!dir.exists(folder)) stop(paste("Folder does not exist:", folder))

  if (is.null(retain)){
    retain <- c("^tracks\\.txt$",                #matches tracks.txt
                #"^slam_trajectory.*\\.txt$",     #matches slam_trajectory... .txt
                "^Accuracy_XY\\.jpg$",           #matches two specific jpgs from the slamreport
                #"^PointCloudAccuracy\\.jpg",     #...
                "^[0-9]{8}.*\\.las$",            #matches the pointcloud.las (another temporary las is excluded by searching for number beginning)
                "SlamReport\\.html$")            #matches the SlamReport
  }

  # Get all files recursively
  files <- list.files(folder, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

  # Identify files to keep
  keep <- c()
  for (f in files) {
    for (p in retain) {
      if (grepl(p, basename(f), perl = TRUE)) {
        keep <- c(keep, f)
        break
      }
    }
  }

  # Move kept files to folder root (flatten structure)
  for (f in keep) {
    dest <- file.path(folder, basename(f))
    if (f != dest) {  # Only move if not already in root
      if (!file.exists(dest)) {
        if (!file.rename(f, dest)) {
          file.copy(f, dest, overwrite = TRUE)
          file.remove(f)
        }
      }
    }
  }

  # Delete files that don't match patterns
  for (f in files) {
    if (!(f %in% keep)) {
      file.remove(f)
    }
  }

  # Delete empty subdirectories (deepest first) - EXCLUDE the folder itself
  dirs <- list.dirs(folder, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[dirs != folder]  # Remove the folder itself from deletion list

  depths <- sapply(dirs, function(d) length(strsplit(d, "/|\\\\")[[1]]))
  dirs <- dirs[order(depths, decreasing = TRUE)]

  for (d in dirs) {
    unlink(d, recursive = TRUE)
  }

  invisible(NULL)
}

#' Retain from a CoPre2 project folder only the relevant files
#'
#' This function performs a standardized filtering of files and subdirectories
#' of the CoPre2 project folder, deletes irrelevant ones and collapses the structure.
#'
#' @param folder   Character. Path to the CoPre2 folder (something like "C:/.../CHC...")
#' @param retain   Character (or Character vector). Regex patterns of relevant
#'                 files to keep. Defaults to
#'                 c("^tracks\\.txt$",
#'                 "^Accuracy_XY\\.jpg$",
#'                 "^[0-9]{8}.*\\.las$",
#'                 "SlamReport\\.html$")
#'
#' @return Returns NULL The actual output is the manipulation of the files and
#'         directories on disk.
#'
#' @export
#' @examples
#' \dontrun{
#' clean_copre_folder("C:/MLS/costa_rica/CHC_20261021111716")
#' }
clean_copre_folder <- function(parent_dir, retain=NULL){

  # Sanity checks
  if (is.null(parent_dir)) stop('Must provide a folder name (something like "C:/.../CHC...")')
  if (!dir.exists(parent_dir)) stop(paste("Parent directory does not exist:", parent_dir))

  # Get scan folders
  scan_folders <- list.files(parent_dir, include.dirs = TRUE, full.names = FALSE, recursive = TRUE)
  scan_folders <- scan_folders[dir.exists(file.path(parent_dir, scan_folders))]
  scan_folders <- scan_folders[dirname(scan_folders) == "Results" &
                                 grepl("^\\d{4}-\\d{2}-\\d{2}-\\d{6}", basename(scan_folders))]

  # Ask to make sure (defaults to TRUE if not interactive)
  if (interactive()) {
    m <- paste0("These are the scan folders in the directory:\n", paste(scan_folders, collapse=", "), "\nAre you sure to proceed?")
    ans <- utils::askYesNo(m)
  } else ans <- TRUE

  # Clean
  if (isTRUE(ans)) {
    # Loop through the scan folders and clean them ####
    for (folder in scan_folders) {
      folder_path <- file.path(parent_dir, folder)
      # Clean the scan folder of unnecessary files ####
      clean_copre_scan(folder_path, retain)
    }
    # Delete further unnecessary files in the project folder ###################
    results <- file.path(parent_dir, "Results")
    stopifnot(dir.exists(results))

    ## 1. Delete all .cpr / .db / .crd files (recursive, case-insensitive) ####
    junk <- list.files(parent_dir, pattern = "\\.(cpr|db|crd)$",
                         recursive = TRUE, full.names = TRUE,
                         ignore.case = TRUE)
    if (length(junk)) unlink(junk, force = TRUE)

    ## 2. Delete log/ and Vector  ####
    unlink(file.path(parent_dir, c("log", "Vector")), recursive = TRUE, force = TRUE)

    ## 3. Move every child of Results/ up one level ####
    items <- list.files(results, full.names = TRUE)
    ok    <- file.rename(items, file.path(parent, basename(items)))
    if (!all(ok)) stop("Failed to move: ", paste(basename(items)[!ok], collapse = ", "))

    ## 4. Remove the now-empty Results/ ####
    unlink(results, recursive = TRUE)
    if (dir.exists(results)) stop("Results/ still exists — check it is empty.")

  } else warning("Cleaning aborted.")
}
