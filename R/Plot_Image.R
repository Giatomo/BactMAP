## collection of raw image plotting functions for BactMAP

# upload your image stack (one color)
#' @export
extr_OriginalStack <- function(picloc) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    inp <- readline("Package 'tiff' needed for this function to work. Press 'y' to install it, or any other key to cancel.")
    if (inp %!in% c("y", "Y")) {
      stop("Canceled")
    }
    utils::install.packages("tiff")
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    inp <- readline("Package 'raster' needed for this function to work. Press 'y' to install it, or any other key to cancel.")
    if (inp %!in% c("y", "Y")) {
      stop("Canceled")
    }
  }
  suppressWarnings(im <- tiff::readTIFF(picloc, all = T)) # if you want the best resolution, it needs to be a .tiff file
  im <- lapply(im, function(x) raster::raster(x))
  imdatframe <- lapply(im, function(x) as.data.frame(methods::as(x, "SpatialPixelsDataFrame"))) # get values
  imdatframe <- lapply(imdatframe, function(x) changecols(x))
  nx <- length(unique(imdatframe[[1]]$x))
  ny <- length(unique(imdatframe[[1]]$y))
  imdatframe <- lapply(imdatframe, function(x) changeres(x, nx, ny))
  return(imdatframe)
}


#' @export
extr_OriginalCells <- function(imdatframe, mesh, surroundings = FALSE, turnCell = TRUE) {
  if ("area" %in% colnames(mesh)) {
    mesh <- mesh[mesh$area > 2, ]
  }

  allcellslist <- lapply(unique(mesh$frame), function(x) pipperframe(imdatframe, mesh, x, surroundings))
  allcellsframe <- do.call(rbind, allcellslist)
  allcellsframe$cell <- allcellsframe$pip
  allcellsframe$pip <- NULL
  if (turnCell == TRUE) {
    outlist <- meshTurn(mesh, rawdatafile = allcellsframe)
  } else {
    outlist <- allcellsframe
  }
  return(outlist)
}


#' @export
plotCellsTime <- function(celdat,
                          updown = T,
                          movie = F,
                          viridisoption = "magma",
                          cellN,
                          minf,
                          maxf,
                          outlines = FALSE,
                          meshdata # ,
                          # overlay=FALSE
) {
  if (movie == TRUE) {
    if (!requireNamespace("gganimate", quietly = TRUE)) {
      inp <- readline("Package 'gganimate' needed to make an animation. Press 'y' to install it, or any other key to cancel.")
      if (inp %in% c("y", "Y")) {
        utils::install.packages("gganimate")
      } else {
        stop("Canceled")
      }
    }
  }
  if (missing(minf)) {
    minf <- min(celdat$frame)
  }
  if (missing(maxf)) {
    maxf <- max(celdat$frame)
  }

  if (outlines == TRUE & !missing(meshdata)) {
    if ("numpoint" %in% colnames(meshdata)) {
      meshdata$num <- meshdata$numpoint
    }
    meshdata <- meshdata[, c("frame", "cell", "X_rot", "Y_rot", "num")]
  }
  if (outlines == TRUE & missing(meshdata)) {
    stop("Outlines is set to TRUE but no meshdata was found. Please specify your mesh dataframe.")
  }
  if (outlines == FALSE & !missing(meshdata)) {
    warning("Meshdata was specified, but outlines are set to FALSE. No outlines will be drawn.")
  }
  # when no cell number is indicated:return a list of plots/movie objects
  if (missing(cellN)) {
    plotout <- lapply(unique(celdat$cell), function(x) {
      plotcellsframelist(celdat[celdat$cell == x, ],
        maxframes = maxf, minframes = minf, updown, movie, viridisoption, outlines, meshdata # , overlay
      ) + ggplot2::ggtitle(x)
    })
  }
  # and just plot one plot if cellN exists
  if (missing(cellN) != T) {
    if (length(cellN) == 1) {
      plotout <- plotcellsframelist(
        celdat[celdat$cell == cellN, ], maxf, minf, updown, movie, viridisoption, outlines, meshdata # , overlay
      ) + ggplot2::ggtitle(cellN)
    }
    if (length(cellN) > 1) {
      plotout <- lapply(cellN, function(x) {
        plotcellsframelist(
          celdat[celdat$cell == x, ], maxf, minf, updown, movie, viridisoption, outlines, meshdata # , overlay
        ) + ggplot2::ggtitle(x)
      })
    }
  }
  return(plotout)
}

################################### |(*-*)/#############################################################
changecols <- function(x) {
  colnames(x) <- c("values", "x", "y")
  return(x)
}
changeres <- function(dat, nx, ny) {
  dat$x <- dat$x * nx - 0.5
  dat$y <- max(dat$y) - dat$y
  dat$y <- dat$y * ny + 1
  return(dat)
}



pinping <- function(dat, mesh, x, surroundings = FALSE) {
  message(paste("Cell", x))
  mesh <- mesh[mesh$cell == x, ]

  if (surroundings == FALSE) {
    minmeshx <- min(mesh$X) - 2
    minmeshy <- min(mesh$Y) - 2
    maxmeshx <- max(mesh$X) + 2
    maxmeshy <- max(mesh$Y) + 2
    dat <- dat[dat$x > minmeshx & dat$y > minmeshy & dat$x < maxmeshx & dat$y < maxmeshy, ]
    # p <- SDMTools::pnt.in.poly(dat[,c("x","y")], mesh[mesh$cell==x,][,c("X","Y")])

    p <- suppressWarnings(sp::point.in.polygon(dat[, "x"], dat[, "y"], mesh$X[mesh$cell == x], mesh$Y[mesh$cell == x])) # find spot/object coordinates inside cell
    p <- data.frame("x" = dat$x, "y" = dat$y, "pip" = p)
    # if pip == 1, the point is inside the polygon. if p==0, it is not.
    # I replace the pips which are 1 with the cell number x
    p$pip[p$pip != 0] <- x
    datje <- merge(dat, p[p$pip != 0, ])
  }
  if (surroundings == TRUE) {
    mesh <- unique(mesh[, c("Xmid", "Ymid", "max.length", "max.width", "angle")])
    w <- 2
    ybox1 <- (0.5 * mesh$max.length + w) * sin(pi - mesh$angle) - (0.5 * mesh$max.width + w) * cos(pi - mesh$angle) + mesh$Ymid
    xbox1 <- (0.5 * mesh$max.length + w) * cos(pi - mesh$angle) + (0.5 * mesh$max.width + w) * sin(pi - mesh$angle) + mesh$Xmid
    xbox2 <- (0.5 * mesh$max.length + w) * cos(pi - mesh$angle) - (0.5 * mesh$max.width + w) * cos(pi - mesh$angle) + mesh$Xmid
    ybox2 <- (0.5 * mesh$max.length + w) * sin(pi - mesh$angle) + (0.5 * mesh$max.width + w) * sin(pi - mesh$angle) + mesh$Ymid
    xbox4 <- mesh$Xmid - (xbox2 - mesh$Xmid)
    ybox4 <- mesh$Ymid - (ybox2 - mesh$Ymid)
    xbox3 <- mesh$Xmid - (xbox1 - mesh$Xmid)
    ybox3 <- mesh$Ymid - (ybox1 - mesh$Ymid)
    box <- data.frame("x" = c(xbox1, xbox2, xbox3, xbox4), "y" = c(ybox1, ybox2, ybox3, ybox4))
    dat <- dat[dat$x > min(box$x) & dat$y > min(box$y) & dat$x < max(box$x) & dat$y < max(box$y), ]
    # p <- SDMTools::pnt.in.poly(dat[,c("x","y")], box)
    p <- suppressWarnings(sp::point.in.polygon(dat[, "x"], dat[, "y"], box$x, box$y))
    p <- data.frame("x" = dat$x, "y" = dat$y, "pip" = p)
    # add cell number (as pip so it matches the situation with the cell outlines)
    p$pip[p$pip != 0] <- x
    datje <- merge(dat, p[p$pip != 0, ])
  }
  return(datje)
}

pipperframe <- function(dat, mesh, y, surroundings = FALSE) {
  mesh <- mesh[mesh$frame == y, ]
  dat <- dat[[y]]
  message(paste("Finding & saving the raw data per cell for frame", y))
  datjeslist <- lapply(unique(mesh$cell), function(x) pinping(dat, mesh, x, surroundings))
  datjesframe <- do.call(rbind, datjeslist)
  datjesframe$frame <- y
  return(datjesframe)
}

################## Plotting cells in a tower/row/movie per cell, per frame

plotcellsframelist <- function(TRframe, maxframes, minframes, updown = F, movie = F, viridisoption = "magma", outlines = FALSE, meshdata # , overlay=FALSE
) {
  nframes <- length(TRframe$frame)
  if (missing(minframes)) {
    minframes <- min(TRframe$frame)
  }
  if (missing(maxframes)) {
    maxframes <- max(TRframe$frame)
  }
  if (nframes > 0) {
    minTRframe <- min(TRframe$frame)
    maxTRframe <- max(TRframe$frame)
    if (minTRframe > minframes) {
      minframes <- minTRframe
    }
    if (maxframes > maxTRframe) {
      maxframes <- maxTRframe
    }

    TRframe <- TRframe[TRframe$frame >= minframes & TRframe$frame <= maxframes, ]

    if (!missing(meshdata) & outlines == TRUE) {
      meshdata <- meshdata[meshdata$frame >= minframes & meshdata$frame <= maxframes, ]
    }

    p <- ggplot2::ggplot(TRframe) +
      ggplot2::geom_polygon(ggplot2::aes_string(x = "xt", y = "yt", fill = "values", group = "pointN"), color = NA) +
      ggimage::theme_transparent() +
      ggplot2::coord_fixed() +
      ggplot2::scale_fill_viridis_c(option = viridisoption) +
      ggplot2::xlim(c(min(TRframe$xt, na.rm = TRUE), max(TRframe$xt, na.rm = TRUE))) +
      ggplot2::ylim(c(min(TRframe$yt, na.rm = TRUE), max(TRframe$yt, na.rm = TRUE))) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        strip.text = ggplot2::element_blank(),
        legend.position = "none",
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0, 0, 0, 0), "mm"),
        panel.spacing = ggplot2::unit(0, "mm"),
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank()
      )

    # if(overlay==TRUE){
    # p <- p + ggplot2::geom_path(data=meshdata, ggplot2::aes(x=X_rot,y=Y_rot), color="white")
    # }

    if (updown == T & movie == F) {
      p <- p + ggplot2::facet_grid(frame ~ .)
    }

    if (updown == F & movie == F) {
      p <- p + ggplot2::facet_grid(~frame)
    }

    if (movie == T) {
      frame <- NULL
      p <- p + gganimate::transition_manual(frame)
    }
    return(p)
  }

  if (nframes == 0) {
    return(NA)
  }
}

### plot raw image

#' @export
plotRaw <- function(tiffdata,
                    meshdata,
                    spotdata,
                    frameN = 1,
                    xrange,
                    yrange,
                    viridisoption = "inferno",
                    meshcolor = "white",
                    spotcolor = "yellow",
                    valuerange,
                    legend = FALSE) {
  if (missing(valuerange) != T & missing(tiffdata) != T) {
    tiffdata[[frameN]] <- tiffdata[[frameN]][tiffdata[[frameN]]$value > valuerange[1] & tiffdata[[frameN]]$value < valuerange[2], ]
  }
  if (legend == TRUE) {
    lpos <- "right"
  }
  if (legend == FALSE) {
    lpos <- "none"
  }
  if (missing(tiffdata) != T) {
    plotcells <- ggplot2::ggplot(tiffdata[[frameN]]) + # plot raw image
      ggplot2::geom_raster(ggplot2::aes_string(x = "x", y = "y", fill = "values")) + # use geom_raster to remake image out of dataframe
      ggplot2::theme_classic() + # simple theme, no backgrounds
      ggplot2::scale_fill_viridis_c(option = viridisoption) + # well-working color scheme for gradient values
      ggplot2::theme(legend.position = lpos) # remove legend for easy viewing
  }
  if (missing(tiffdata) == T) {
    plotcells <- ggplot2::ggplot() +
      ggplot2::theme_dark()
  }
  # add x and/or y range + fixed coordinates when indicated:
  if (missing(xrange) != T & missing(yrange) != T) {
    plotcells <- plotcells + ggplot2::coord_fixed(xlim = xrange, ylim = yrange) # sub-set of the image frame to zoom in
  }
  if (missing(xrange) != T & missing(yrange) == T) {
    plotcells <- plotcells + ggplot2::coord_fixed(xlim = xrange)
  }
  if (missing(xrange) == T & missing(yrange) != T) {
    plotcells <- plotcells + ggplot2::coord_fixed(ylim = yrange)
  }
  if (missing(xrange) == T & missing(yrange) == T) {
    plotcells <- plotcells + ggplot2::coord_fixed()
  }


  # add mesh data when given:
  if (missing(meshdata) != T) {
    meshdata <- meshdata[order(meshdata$frame, meshdata$cell, meshdata$num), ]
    plotcells <- plotcells + # plot made above
      ggplot2::geom_path(data = meshdata[meshdata$frame == frameN, ], ggplot2::aes_string(x = "X", y = "Y", group = "cell"), color = meshcolor) # add outline of cells, only frame one, white color
  }
  if (missing(spotdata) != T) {
    plotcells <- plotcells +
      ggplot2::geom_point(data = spotdata[spotdata$frame == frameN, ], ggplot2::aes_string(x = "x", y = "y"), shape = 1, color = spotcolor) # add yellow empty dots of our spot localizations on top
  }
  return(plotcells)
}



st_box <- function(left, right, top, bot, center = NULL, width = NULL, height = NULL) {
  if (!missing(center)) {
    if (missing(width) | missing(height)) {
      stop("you need to provide 'width' and 'height' arguments when using 'center' argument")
    }
    if (missing(left) | missing(right)) {
      stop("you need to provide 'width' and 'height' arguments when using 'center' argument")
    }
    left <- center - width/2
    right <- center + width/2
    bot <- center - height/2
    top <- center + height/2
  }
  if (xor(missing(top) | missing(right), (missing(width) | missing(heigth)))) {
    stop("you need to provide 'width' and 'height' OR 'top' arguments when using 'left' and bot argument")
    }
  if (missing(top) | missing(right)) {
    right <- left + width
    top <- bot + heigth
  }
  
  points <- matrix(c(right, top, left, top, left, bot, right, bot, right, top), ncol = 2, byrow = TRUE)
  return(st_polygon(list(points)))
}





######## plot time series (or by cell length) raw data kymographs
#' @export
TurnedCell4 |> 
  select(c(values, frame, cell, X_rot, Y_rot, max.length, max.width, pointN)) |> 
  distinct() |>
  group_by(frame, cell) |> 
  mutate(
    group = cut(X_rot, breaks = bins, labels = 1:bins),
#    group_width  = cut(Y_rot, breaks = bins, labels = 1:bins)
    ) |>
  group_by(frame, cell, group) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
  group_by(frame, cell) |>
  mutate(
    right = frame + 0.5, #a
    left = frame - 0.5, #b
    bot = X_rot - (max(X_rot) / bins) * 1.15, #c
    top = X_rot + (max(X_rot) / bins) * 1.15) |> #d
  pivot_longer(c(right, left), names_to = "frame_location", values_to = "x_coords") |>
  pivot_longer(c(bot, top), names_to = "variable", values_to = "y_coords") |>
  mutate(
    frame_location = if_else(frame_location == "left" & variable == "bot", "bot", frame_location),
    frame_location = if_else(frame_location == "right" & variable == "bot", "top", frame_location)
    ) |>
  select(-variable) |>
  ungroup() |>
  arrange(max.length, frame, cell, group, frame_location) |>
#  group_by(cell, frame, max.width) |>
#  mutate(id = cur_group_id()) |>
#  group_by(group, frame_location, frame) |> 
#  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
  group_by(cell, frame, X_rot, pointN, values) |>
  mutate(id = cur_group_id()) -> original
  
  # using sf and polygons (plot a bit nicer imo)
  
  TurnedCell4 |> 
  select(c(values, frame, cell, X_rot, Y_rot, max.length, max.width, pointN)) |> 
  distinct() |>
  group_by(frame, cell) |> 
  mutate(
    group = cut(X_rot, breaks = bins, labels = 1:bins),
#    group_width  = cut(Y_rot, breaks = bins, labels = 1:bins)
#    ) |>
  group_by(frame, cell, group) |>
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
  group_by(frame, cell) |>
  mutate(
    right = frame + 0.5, 
    left = frame - 0.5, 
    bot = X_rot - (max(X_rot) / bins) * 1.15, 
    top = X_rot + (max(X_rot) / bins) * 1.15,
    ) |>
  rowwise() |>
  arrange(max.length, frame, cell, group) |>
  mutate(polygon = list(st_box(left, right, top, bot))) |> st_sf() -> with_sf
  

  with_sf |> ggplot() + geom_sf(aes(fill = values, color = values))
  original |> ggplot() + geom_polygon(aes(x = x_coords, y = y_coords, group = id, fill = values, color = values))
  
prepForKymo <- function(turnedCells, dimension = "length", bins = 25, sizeAV = FALSE) {
  turnedCells <- turnedCells |> 
    select(c(values, frame, cell, X_rot, Y_rot, max.length, max.width, pointN)) |> 
    distinct() 
    
  originalCells <- list()
  for (n in unique(turnedCells$frame)) {
    U <- turnedCells[turnedCells$frame == n, ]
    if (dimension == "length") {
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$X_rot[U$cell == x], breaks = bins, labels = c(1:bins))))
    }
    if (dimension == "width") {
      U$group <- unlist(lapply(unique(U$cell), function(x) cut(U$Y_rot[U$cell == x], breaks = bins, labels = c(1:bins))))
    }

    Umeans <- suppressWarnings(lapply(unique(U$cell), function(x) stats::aggregate(U[U$cell == x, ], by = list(U$group[U$cell == x]), FUN = mean, na.rm = T)))

    Umeansall <- do.call("rbind", Umeans)


    if (sizeAV == TRUE) {
      Umeansall$a <- Umeansall$frame + 0.5
      Umeansall$b <- Umeansall$frame - 0.5
      if (dimension == "length") {
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell == x] + (max(Umeansall$X_rot[Umeansall$cell == x]) / bins) * 1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$X_rot[Umeansall$cell == x] - (max(Umeansall$X_rot[Umeansall$cell == x]) / bins) * 1.15))
      }
      if (dimension == "width") {
        Umeansall$d <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell == x] + (max(Umeansall$Y_rot[Umeansall$cell == x]) / bins) * 1.15))
        Umeansall$c <- unlist(lapply(unique(Umeansall$cell), function(x) Umeansall$Y_rot[Umeansall$cell == x] - (max(Umeansall$Y_rot[Umeansall$cell == x]) / bins) * 1.15))
      }


      Umeansall <- reshape2::melt(Umeansall, measure.vars = c("a", "b"), value.name = "x_coords")
      Umeansall$frameh <- Umeansall$variable
      Umeansall$variable <- NULL
      Umeansall <- reshape2::melt(Umeansall, measure.vars = c("d", "c"), value.name = "y_coords")
      Umeansall$frameh <- as.character(Umeansall$frameh)
      Umeansall$frameh[Umeansall$frameh == "b" & Umeansall$variable == "c"] <- "c"
      Umeansall$frameh[Umeansall$frameh == "a" & Umeansall$variable == "c"] <- "d"
      Umeansall$variable <- NULL
      Umeansall <- Umeansall[order(Umeansall$frame, Umeansall$cell, Umeansall$Group.1, Umeansall$frameh), ]
    }

    Umeansall$group <- Umeansall$Group.1
    Umeansall$Group.1 <- NULL

    originalCells[[n]] <- Umeansall
  }

  originalCells <- do.call("rbind", originalCells)

  originalCells <- originalCells[order(originalCells$max.length), ]
  cellnumdatframe <- data.frame(cellnum.length = 1:length(unique(originalCells$max.length)), max.length = unique(originalCells$max.length))
  originalCells <- merge(originalCells, cellnumdatframe)

  originalCells <- originalCells[order(originalCells$max.width), ]
  cellnumdatframe <- data.frame(cellnum.width = 1:length(unique(originalCells$max.width)), max.width = unique(originalCells$max.width))
  originalCells <- merge(originalCells, cellnumdatframe)

  return(originalCells)
}


############## plot#########
#' @export

bactKymo <- function(originalCells, 
                     timeD = FALSE, dimension = "length",
                      bins = 25, sizeAV = FALSE,
                      cells = "all", prep = TRUE,
                      percDiv = FALSE, cutoff_demograph = 0.975, 
                      mag, legend = TRUE) {
  
  measure <- "pixels"
  pos <- if_else(legend, "right", "none")
  cells <- ifelse(
    cells[1] == "all",
    unique(originalCells$cell[order(originalCells$cell)]),
    cells)
  
  if (!is.numeric(cells)) {
      stop("'cells' must be either be a character string 'all', a number or a vector of numbers c(x, x,)")
    }
  
  if (percDiv) {
    if ("percentage_binned" %in% colnames(originalCells)) {
      groupP <- unique(originalCells[, c("cell", "frame", "percentage_binned")])
    }
    else {
      groupP <- perc_Division(unique(originalCells[, c("cell", "frame", "max.length")]), av = FALSE, plotgrowth = FALSE)$timelapse
      groupP <- unique(groupP[, c("cell", "frame", "percentage_binned")])
    }
  }

  if (prep == TRUE) {
    originalCells <- prepForKymo(originalCells, dimension = dimension, bins = bins, sizeAV = sizeAV)
    if (percDiv == TRUE) {
      originalCells <- merge(originalCells, groupP)
    }
  }

  if (percDiv) {
    originalCells$frame <- as.numeric(originalCells$percentage_binned)
    originalCells$cell <- 1
    
    
    if (sizeAV) {
      originalCells <- stats::aggregate(originalCells[, colnames(originalCells)[colnames(originalCells) != "group" & colnames(originalCells) != "frameh" & colnames(originalCells) != "frame" & colnames(originalCells) != "percentage_binned"]],
        by = list(group = originalCells$group, frameh = originalCells$frameh, frame = originalCells$frame), FUN = mean
      )

      originalCells$x_coords[originalCells$frameh == "a" | originalCells$frameh == "d"] <- originalCells$frame[originalCells$frameh == "a" | originalCells$frameh == "d"] + 0.5
      originalCells$x_coords[originalCells$frameh == "b" | originalCells$frameh == "c"] <- originalCells$frame[originalCells$frameh == "b" | originalCells$frameh == "c"] - 0.5
    }

    else {
      originalCells <- stats::aggregate(originalCells[, colnames(originalCells)[colnames(originalCells) != "group" & colnames(originalCells) != "frame" & colnames(originalCells) != "percentage_binned"]],
        by = list(group = originalCells$group, frame = originalCells$frame), FUN = mean
      )
    }
    
  }
  
  if (!missing(mag)) {
        originalCells$y_coords <- originalCells$y_coords * unlist(get(magnificationList, envir = magEnv)[mag])
        originalCells$x_coords <- originalCells$x_coords * unlist(get(magnificationList, envir=magEnv)[mag])
        measure <- "micron"
      }
  
  if (percDiv) {
    if (sizeAV) {
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$values, sep = "_")
      originalCells$x_coords <- originalCells$x_coords * 10

      plot1 <- ggplot2::ggplot(originalCells) +
        ggplot2::geom_polygon(ggplot2::aes_string(x = "x_coords", y = "y_coords", fill = "values", group = "grouping")) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Percentage of division") +
        ggplot2::ylab(paste("length by cell ", dimension, "( in", measure, ")", sep = "")) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme(legend.position = pos)
        
      return (plot1)
    }
    
    plot1 <- ggplot2::ggplot(originalCells) +
      ggplot2::geom_raster(ggplot2::aes_string(x = "frame", y = "group", fill = "values")) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("Percentage of division") +
      ggplot2::ylab(paste("length (by cell ", dimension, ")", sep = "")) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme(legend.position = pos)
      
      return (plot1)
    }

  if (timeD) {
    if (sizeAV) {
      
      originalCells$grouping <- paste(originalCells$X_rot, originalCells$pointN, originalCells$values, sep = "_")
      plot1 <- lapply(cells, function(x) {
        ggplot2::ggplot(originalCells[originalCells$cell == x, ]) +
          ggplot2::geom_polygon(ggplot2::aes_string(x = "x_coords", y = "y_coords", fill = "values", group = "grouping")) +
          ggplot2::theme_minimal() +
          ggplot2::xlab("Time (frames)") +
          ggplot2::ylab(paste(dimension, " (in ", measure, ")", sep = "")) +
          ggplot2::scale_fill_viridis_c() +
          ggplot2::ggtitle(paste("Cell", x, sep = " ")) +
          ggplot2::theme(legend.position = pos)
      })
      
      if (length(cells) == 1) {
          plot1 <- plot1[[1]]
          return (plot1)
      }
      
      names(plot1) <- paste("cell", cells, sep = "")
      return (plot1)

  }
    plot1 <- lapply(cells, function(x) {
      ggplot2::ggplot(originalCells[originalCells$cell == x, ]) +
        ggplot2::geom_raster(ggplot2::aes_string(x = "frame", y = "group", fill = "values")) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Time (frames)") +
        ggplot2::ylab(paste("bin (by cell ", dimension, ")", sep = "")) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::ggtitle(paste("Cell", x, sep = " ")) +
        ggplot2::theme(legend.position = pos)
    })
    
    if (length(cells) == 1) {
      plot1 <- plot1[[1]]
      return(plot1)
    }
    
    names(plot1) <- paste("cell", cells, sep = "")
    return (plot1)
  }

  # 97.5 % cutoff for outlying superbright stuff
  if (length(cells) <= 1) {
    stop("'cells' is of length 1, while there is no time dimension.
          \n if you want to plot a single cell over time, set 'timeD' to 'TRUE'
          \n if you want to plot all cells as a demograph, put 'cells' to 'all'.
          \n if you want to plot a specific group of cells as a demograph, put 'cells' to a vector identifying the cell numbers")
  }

  if (!is.numeric(cells)) {
    stop("'cells' is neither numeric nor 'all', please set 'cells' to a numeric vector or 'all'")
  }
  
  originalCells <- originalCells[originalCells$cell %in% cells, ]
  originalCells <- originalCells[originalCells$values < stats::quantile(originalCells$values, cutoff_demograph), ]

  cellnum.dimension <- switch(dimension,
         "length" = "cellnum.length",
         "width" = "cellnum.width"
         )

  if (!sizeAV ) {
    plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x = cellnum.dimension, y = "group", fill = "values")) +
      ggplot2::geom_raster() +
      ggplot2::coord_fixed(ratio = 20) +
      ggplot2::theme_minimal() +
      ggplot2::xlab(paste("n(th) cell ordered by cell", dimension)) +
      ggplot2::ylab(paste("bin (by cell", dimension, ")")) +
      ggplot2::scale_fill_viridis_c(name = "Fluorescence\nIntensity") +
      ggplot2::theme(legend.position = pos)
    return (plot1)
    }
  
    originalCells[cellnum.dimension][originalCells$frameh == "a" | originalCells$frameh == "d"] <- originalCells[cellnum.dimension][originalCells$frameh == "a" | originalCells$frameh == "d"] + 0.5
    originalCells[cellnum.dimension][originalCells$frameh == "b" | originalCells$frameh == "c"] <- originalCells[cellnum.dimension][originalCells$frameh == "b" | originalCells$frameh == "c"] - 0.5
    
    originalCells$grouping <- paste(originalCells$X_rot, originalCells$cell, originalCells$frame, originalCells$pointN, sep = "_")
    plot1 <- ggplot2::ggplot(originalCells, ggplot2::aes_string(x = cellnum.dimension, y = "y_coords", group = "grouping", fill = "values")) +
      ggplot2::geom_polygon() +
      ggplot2::theme_minimal() +
      ggplot2::xlab(paste("n(th) cell ordered by cell", dimension) +
      ggplot2::ylab(paste("location on length axis (", measure, ")", sep = "")) +
      ggplot2::scale_fill_viridis_c(name = "Fluorescence\nIntensity") +
      ggplot2::theme(legend.position = pos)
    return(plot1)
    }

                    
GeomKymo <- ggplot2::ggproto(
  "GeomKymo",
  c(ggplot2::Geom, ggplot2::Stat,
  compute_group = function(data, scales) {
    data |> 
      arrange(y) |> 
      group_by(x) |>
      mutate(signal = list(sample(100, y)),
             y_corr = y - y/2)
  }
))

geom_kymo <- ggplot2




      
      
      
      