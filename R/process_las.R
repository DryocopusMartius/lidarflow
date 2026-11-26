###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
#      name:  minimal_processing_lasR                                          #
# timestamp:  2026-10-06T12:46Z                                                #
#    author:  Tillman Reuter (orcid.org/0009-0002-1156-9983)                   #
#   license:  GPLv3 (see LICENSE.md)                                           #
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%

# READ ME ######################################################################
'
 PURPOSE:
   This script provides a robust, documented function for processing LiDAR
   point clouds using the lasR package. It serves as a template for
   implementing a 2-stage TIN (Triangulated Irregular Network) approach
   for ground classification and normalization.

 DEPENDENCIES:
   - lasR   : Core processing engine
   - lidR   : Used for specific pre-clipping tasks (trajectory hull)
   - sf     : Required for spatial operations during trajectory clipping
   - parallel : Used for multi-core processing detection

 USAGE:
   Source this file in your R session, then call `process_las()`.
'
par(mfrow=c(1,1))
Sys.setenv(LANG = "en")
#Sys.setenv(LANG = "ger")

###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
# LOAD AND INSTALL PACKAGES ####################################################
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%

# Check and load required libraries. If missing, installation is suggested.
required_packages <- c("lasR", "lidR", "sf")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Package '%s' is not installed. Please install it manually.", pkg))
  }
}

###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%
# FUNCTIONS ####################################################################
###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%###%


#' Process LiDAR Point Clouds using lasR
#'
#' This function performs a standardized processing workflow on LAS/LAZ files.
#' It supports clipping, denoising, subsampling, ground classification,
#' normalization, and derivative product generation (DTM, height slices).
#' It utilizes the \code{lasR} package for the main pipeline and \code{lidR}
#' for specific pre-processing tasks like trajectory-based clipping.
#'
#' @param input_file   Character. Path to the input .las or .laz file.
#' @param output_dir   Character. Path to the directory where output files
#'                     will be saved. If \code{NULL}, defaults to the directory
#'                     of the input file.
#' @param drop_color   Logical. If \code{TRUE}, RGB color information is
#'                     discarded to save compute. Defaults to \code{FALSE}.
#' @param clip_to      Character. Defines the clipping strategy. Options:
#'                     \itemize{
#'                       \item \code{"none"}: No clipping (default).
#'                       \item \code{"circle"}: Clip to a circular area.
#'                       \item \code{"traj"}: Clip to the buffered convex hull
#'                       of a trajectory file.
#'                     }
#' @param x_center     Numeric. X coordinate of the circle center (required if
#'                     \code{clip_to = "circle"}).
#' @param y_center     Numeric. Y coordinate of the circle center (required if
#'                     \code{clip_to = "circle"}).
#' @param radius       Numeric. Radius of the clipping circle in meters
#'                     (required if \code{clip_to = "circle"}).
#' @param traj_path    Character. Path to a CSV file containing trajectory
#'                     points (required if \code{clip_to = "traj"}).
#'                     Expected columns: ID, X, Y.
#' @param buffer       Numeric. Buffer distance in meters to apply around the
#'                     trajectory hull.
#' @param denoise      Logical. If \code{TRUE}, applies Statistical Outlier
#'                     Removal (SOR) to classify and remove noise.
#' @param sor_k        Integer. Number of neighbors for SOR filter.
#'                     Defaults to package standard if \code{NULL}.
#' @param sor_m        Numeric. Standard deviation multiplier for SOR filter.
#' @param subsample    Logical. If \code{TRUE}, applies Poisson disk subsampling.
#' @param distance     Numeric. Minimum 3D distance between output points
#'                     (meters) for subsampling.
#' @param normalize    Logical. If \code{TRUE}, performs ground classification
#'                     and normalizes heights (Z relative to ground).
#' @param start_res    Numeric. Initial resolution for Progressive TIN Densification (PTD).
#' @param angle        Numeric. Maximum angle (degrees) for PTD iteration.
#' @param output_dtm   Logical. If \code{TRUE}, generates and saves a Digital
#'                     Terrain Model (GeoTIFF).
#' @param output_below1m Logical. If \code{TRUE}, saves a separate point cloud
#'                     containing only points below 1 meter height.
#' @param output_unnormalized Logical. If \code{TRUE}, saves a copy of the
#'                     point cloud before normalization.
#' @param dtm_res      Numeric. Resolution of the output DTM in meters.
#' @param ncores       Integer. Number of CPU cores to use. If \code{NULL},
#'                     automatically detects and uses 90\% of available cores.
#' @param verbose      Logical. If \code{TRUE}, prints progress messages to
#'                     the console.
#'
#' @return Invisibly returns the result of the \code{lasR::exec} call.
#'         Output files are written to disk as side effects.
#'
#' @export
process_las <- function(input_file,
                        output_dir         = NULL,

                        drop_color         = FALSE,

                        clip_to            = "none",
                        x_center           = NULL,
                        y_center           = NULL,
                        radius             = NULL,
                        traj_path          = NULL,
                        buffer             = NULL,

                        denoise            = FALSE,
                        sor_k              = NULL,
                        sor_m              = NULL,

                        subsample          = FALSE,
                        distance           = NULL,

                        normalize          = FALSE,
                        start_res          = 1,
                        angle              = 30,

                        output_dtm         = FALSE,
                        output_below1m     = FALSE,
                        output_unnormalized = FALSE,
                        dtm_res            = 0.5,

                        ncores             = NULL,
                        verbose            = TRUE) {

  # ============================================================================
  # INITIALIZATION AND SANITY CHECKS
  # ============================================================================

  begin_time <- Sys.time()

  # Warn if no processing steps are selected
  if (all(!drop_color, !denoise, !subsample, !normalize, clip_to == "none")) {
    warning("No meaningful processing steps are performed. Check arguments.")
  }

  # Validate input file existence
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  # Validate circle clipping parameters
  if (clip_to == "circle") {
    if (is.null(x_center) | is.null(y_center)) {
      x_center <- y_center <- 0
      message("Circle center is not specified. Defaulting to (0, 0)")
    }
    if (is.null(radius)) stop("Argument 'radius' is required when clip_to = 'circle'.")
  }

  # Validate trajectory clipping parameters
  if (clip_to == "traj") {
    if (is.null(traj_path)) {
      stop("Argument 'traj_path' is required when clip_to = 'traj'.")
    }
    if (!file.exists(traj_path)) {
      stop("Trajectory file not found: ", traj_path)
    }
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for clip_to = 'traj'.")
    }
  }

  # Validate resolution consistency
  if (all(!is.null(distance), !is.null(dtm_res))){
    if (distance>dtm_res) {
      stop(
        "DTM resolution (", dtm_res,
        ") is higher than subsampling distance (", distance,
        "). This is not recommended."
      )
    }
  }

  # ============================================================================
  # OUTPUT PATH CONFIGURATION
  # ============================================================================

  # Determine output directory
  if (is.null(output_dir)){
    message("No output directory specified. Using input file directory.")
    output_dir <- dirname(input_file)
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Construct output file paths
  out_base          <- tools::file_path_sans_ext(basename(input_file))
  file_processed    <- file.path(output_dir, paste0(out_base, "_processed.las"))
  file_unnormalized <- file.path(output_dir, paste0(out_base, "_unnormalized.las"))
  file_normalized   <- file.path(output_dir, paste0(out_base, "_normalized.las"))
  file_cutoff       <- file.path(output_dir, paste0(out_base, "_below_1m.las"))
  file_dtm          <- file.path(output_dir, paste0(out_base, "_dtm.tif"))

  # ============================================================================
  # PRE-PROCESSING: TRAJECTORY CLIPPING (HYBRID lidR APPROACH)
  # ============================================================================

  working_file <- input_file   # defaults to original file

  if (clip_to == "traj") {

    if (verbose) message("Pre-clipping to trajectory hull using lidR and sf...")

    # 1. Read trajectory and create spatial hull
    traj <- read.table(traj_path, sep = ",", header = FALSE)
    traj_xy <- traj[, 2:3]
    traj_xy <- -traj_xy
    names(traj_xy) <- c("X", "Y")

    # Convert to SF object and compute convex hull
    pts_sf <- sf::st_as_sf(traj_xy, coords = c("X", "Y"))
    hull   <- sf::st_convex_hull(sf::st_union(pts_sf))

    # Apply buffer if requested
    if (!is.null(buffer)) hull <- sf::st_buffer(hull, dist = buffer)

    # 2. Read original LAS with lidR and clip to the polygon
    las_full    <- lidR::readLAS(input_file)
    las_clipped <- lidR::clip_roi(las_full, hull)

    if (lidR::npoints(las_clipped) == 0) {
      stop("Trajectory hull contains no points. Check CRS and buffer.")
    }

    # 3. Write clipped cloud to a temporary file for lasR processing
    working_file <- tempfile(fileext = ".las")
    lidR::writeLAS(las_clipped, working_file)

    # Clean memory
    rm(las_full, las_clipped, pts_sf, hull)
    gc()

    if (verbose) message("Clipped cloud written to temporary file: ", working_file)
  }

  # ============================================================================
  # BUILD lasR PIPELINE
  # ============================================================================

  stages <- list()

  # --- 1. Reader ---
  if (clip_to == "circle") {
    if (verbose) message("Adding circle cut (r = ", radius, " m).")
    stages$reader <- lasR::reader_circles(
      xc = x_center,
      yc = y_center,
      r  = radius)
  } else {
    stages$reader <- lasR::reader()
  }

  # --- 2. Dropping color ---  if (drop_color){
  if (verbose){
    message("Dropping color attributes (RGB)")
    stages$drop_color <- lasR::remove_rgb()
  }

  # --- 3. Subsampling ---
  if (subsample) {
    if (verbose) message("Adding Poisson subsampling (", distance, " m).")
    stages$subsample <- lasR::sampling_poisson(distance = distance)
  }

  # --- 4. Denoising ---
  if (denoise) {
    if (verbose) message("Adding noise filter (k=",sor_k,", m=",sor_m,").")
    stages$classify_noise <- lasR::classify_with_sor(k = sor_k,
                                                     m = sor_m,
                                                     class = 18L)
    stages$delete_noise   <- lasR::delete_noise()
  }

  # --- 5. Ground Classification & Normalization ---
  if (normalize) {

    # Optional: Save unnormalized state
    if (output_unnormalized){
      if (verbose)            message("Writing cloud without normalization to: ", file_unnormalized)
      stages$write_unnorm <-  lasR::write_las(file_unnormalized)
    }
    if (verbose)              message("Terrain normalizing (PTD)")

    # Additional noise filter specifically for PTD stability
    stages$noise_for_ptd <-   lasR::classify_with_sor(k=20,m=0, class = 18L)

    # Progressive TIN Densification
    stages$classify_ground <- lasR::classify_with_ptd(res=start_res, angle = angle, filter = lasR::drop_class(18L))

    # Triangulate ground points
    stages$tri <- tri <-      lasR::triangulate(filter = lasR::keep_ground())

    # Optional: Generate DTM Raster
    if (output_dtm){
      if (verbose)            message("Writing DTM to: ", file_dtm)
      stages$dtm <-           lasR::rasterize(res = dtm_res, tri, ofile = file_dtm)
    }
    # Normalize heights using the TIN
    stages$transform <-       lasR::transform_with(tri)

    # Cleanup: Remove negative heights (below ground artifacts)
    stages$lower_cutoff <-    lasR::delete_points("Z < -0.1")

    # Save Normalized Cloud
    if (verbose)              message("Writing normalized cloud to: ", file_normalized)
    stages$write_norm <-      lasR::write_las(file_normalized)

    # Optional: Vertical Cutoff (Vegetation/Understory)
    if (output_below1m){
      if (verbose)            message("Adding vertical cutoff: removing points with Z > 1 m. Writing below_1m cloud to: ", file_cutoff)
      stages$upper_cutoff <-  lasR::write_las(file_cutoff, filter = "Z <= 1")
    }
  } else {
    # If not normalizing, just write the processed result
    if (verbose)              message("Writing output cloud to: ", file_processed)
    stages$write_unnorm <-    lasR::write_las(file_processed)
  }

  # ============================================================================
  # EXECUTE PIPELINE
  # ============================================================================

  # Determine number of cores
  if (is.null(ncores)){
    available_cores <- parallel::detectCores()
    ncores <- ceiling(0.9*available_cores)
    if (verbose) message("Using ", ncores, " out of ", available_cores, " available cores.")
  }

  # Combine all stages into a single pipeline
  pipeline <- Reduce(`+`, stages)

  # Run the pipeline
  result <- lasR::exec(pipeline, on = working_file, ncores = ncores, progress = TRUE)

  # ============================================================================
  # CLEANUP AND FINALIZATION
  # ============================================================================

  # Remove temporary file if trajectory clipping was used
  if (clip_to == "traj" && file.exists(working_file)) {
    unlink(working_file)
  }

  # Report execution time
  end_time <- Sys.time()
  if (verbose) message("Process took ", format(end_time-begin_time))

  # Return result invisibly
  invisible(result)
}

