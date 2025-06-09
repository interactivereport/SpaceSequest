#This is to fix the CosMx data loading error
library(Seurat)
library(progressr)
LoadNanostring_fix <- function(data.dir, fov, assay = 'Nanostring') {
  data <- ReadNanostring_fix(
    data.dir = data.dir,
    type = c("centroids", "segmentations")
  )
  segs <- CreateSegmentation(data$segmentations)
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = segs
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$pixels,
    assay = assay
  )
  obj <- CreateSeuratObject(counts = data$matrix, assay = assay)
  
  # subset both object and coords based on the cells shared by both
  cells <- intersect(
    Cells(x = coords, boundary = "segmentation"),
    Cells(x = coords, boundary = "centroids")
  )
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[[fov]] <- coords
  return(obj)
}


ReadNanostring_fix <- function(
    data.dir,
    mtx.file = NULL,
    metadata.file = NULL,
    molecules.file = NULL,
    segmentations.file = NULL,
    type = 'centroids',
    mol.type = 'pixels',
    metadata = NULL,
    mols.filter = NA_character_,
    genes.filter = NA_character_,
    fov.filter = NULL,
    subset.counts.matrix = NULL,
    cell.mols.only = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('centroids', 'segmentations'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels'),
    several.ok = TRUE
  )
  if (!is.null(metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c(
        "Area", "fov", "Mean.MembraneStain", "Mean.DAPI", "Mean.G",
        "Mean.Y", "Mean.R", "Max.MembraneStain", "Max.DAPI", "Max.G",
        "Max.Y", "Max.R"
      ),
      several.ok = TRUE
    )
  }
  
  use.dir <- all(vapply(
    X = c(mtx.file, metadata.file, molecules.file),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Nanostring directory ", data.dir)
  }
  # Identify input files
  files <- c(
    matrix = mtx.file %||% '[_a-zA-Z0-9]*_exprMat_file.csv',
    metadata.file = metadata.file %||% '[_a-zA-Z0-9]*_metadata_file.csv',
    molecules.file = molecules.file %||% '[_a-zA-Z0-9]*_tx_file.csv',
    segmentations.file = segmentations.file %||% '[_a-zA-Z0-9]*-polygons.csv'
  )
  
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  
  if (all(is.na(x = files))) {
    stop("Cannot find Nanostring input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['metadata.file']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    md <- data.table::fread(
      file = files[['metadata.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      md <- md[md$fov %in% fov.filter,]
    }
    pprecoord(type = 'finish')
  }
  if (!is.na(x = files[['segmentations.file']])) {
    ppresegs <- progressor()
    ppresegs(
      message = "Preloading cell segmentation vertices",
      class = 'sticky',
      amount = 0
    )
    segs <- data.table::fread(
      file = files[['segmentations.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      segs <- segs[segs$fov %in% fov.filter,]
    }
    ppresegs(type = 'finish')
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules.file']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules.file']],
      sep = ',',
      verbose = FALSE
    )
    
    # filter molecules file by FOVs
    if (!is.null(x = fov.filter)) {
      mx <- mx[mx$fov %in% fov.filter,]
    }
    
    # Molecules outside of a cell have a cell_ID of 0
    if (cell.mols.only) {
      mx <- mx[mx$cell_ID != 0,]
    }
    
    if (!is.na(x = mols.filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", mols.filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = mols.filter, x = mx$target), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules.file']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules.file')]
  }
  files <- files[!is.na(x = files)]
  
  outs <- list("matrix"=NULL, "pixels"=NULL, "centroids"=NULL)
  if (!is.null(metadata)) {
    outs <- append(outs, list("metadata" = NULL))
  }
  if ("segmentations" %in% type) {
    outs <- append(outs, list("segmentations" = NULL))
  }
  
  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'matrix' = {
        ptx <- progressor()
        ptx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        if (!is.null(subset.counts.matrix)) {
          tx <- build.cellcomp.matrix(mols.df=mx, class=subset.counts.matrix)
        } else {
          tx <- data.table::fread(
            file = files[[otype]],
            sep = ',',
            data.table = FALSE,
            verbose = FALSE
          )
          # Combination of Cell ID (for non-zero cell_IDs) and FOV are assumed to be unique. Used to create barcodes / rownames.
          bcs <- paste0(as.character(tx$cell_ID), "_", tx$fov)
          rownames(x = tx) <- bcs
          # remove all rows which represent counts of mols not assigned to a cell for each FOV
          tx <- tx[!tx$cell_ID == 0,]
          # filter fovs from counts matrix
          if (!is.null(x = fov.filter)) {
            tx <- tx[tx$fov %in% fov.filter,]
          }
          # atomx export has an additional column called 'cell' Simply remove it!
          # This casuses Seurat function to fail
          # tx <- subset(tx, select = -c(fov, cell_ID))
          tx <- subset(tx, select = -c(fov, cell, cell_ID))
        }
        
        tx <- as.data.frame(t(x = as.matrix(x = tx)))
        if (!is.na(x = genes.filter)) {
          ptx(
            message = paste("Filtering genes with pattern", genes.filter),
            class = 'sticky',
            amount = 0
          )
          tx <- tx[!grepl(pattern = genes.filter, x = rownames(x = tx)), , drop = FALSE]
        }
        # only keep cells with counts greater than 0
        tx <- tx[, which(colSums(tx) != 0)]
        ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
        
        if ((sum(tx == 0) / length(x = tx)) > ratio) {
          ptx(
            message = 'Converting counts to sparse matrix',
            class = 'sticky',
            amount = 0
          )
          tx <- as.sparse(x = tx)
        }
        
        ptx(type = 'finish')
        
        tx
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Creating centroid coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = md$CenterX_global_px,
          y = md$CenterY_global_px,
          cell = paste0(as.character(md$cell_ID), "_", md$fov),
          stringsAsFactors = FALSE
        )
      },
      'segmentations' = {
        pcents <- progressor()
        pcents(
          message = 'Creating segmentation coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = segs$x_global_px,
          y = segs$y_global_px,
          cell = paste0(as.character(segs$cellID), "_", segs$fov),  # cell_ID column in this file doesn't have an underscore
          stringsAsFactors = FALSE
        )
      },
      'metadata' = {
        pmeta <- progressor()
        pmeta(
          message = 'Loading metadata',
          class = 'sticky',
          amount = 0
        )
        pmeta(type = 'finish')
        df <- md[,metadata]
        df$cell <- paste0(as.character(md$cell_ID), "_", md$fov)
        df
      },
      'pixels' = {
        ppixels <- progressor()
        ppixels(
          message = 'Creating pixel-level molecule coordinates',
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$x_global_px,
          y = mx$y_global_px,
          gene = mx$target,
          stringsAsFactors = FALSE
        )
        ppixels(type = 'finish')
        df
      },
      # 'microns' = {
      #   pmicrons <- progressor()
      #   pmicrons(
      #     message = "Creating micron-level molecule coordinates",
      #     class = 'sticky',
      #     amount = 0
      #   )
      #   df <- data.frame(
      #     x = mx$global_x,
      #     y = mx$global_y,
      #     gene = mx$gene,
      #     stringsAsFactors = FALSE
      #   )
      #   pmicrons(type = 'finish')
      #   df
      # },
      stop("Unknown Nanostring input type: ", outs[[otype]])
    )
  }
  return(outs)
}