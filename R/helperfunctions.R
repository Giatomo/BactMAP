## Renske van Raaphorst
## 7/14/2017

## Helperfunctions: functions necessary to make other functions working properly.
## other package dependencies:

# merge data functions

`%!in%` <- Negate(`%in%`)

#' @export
#'
#' @title Function to get a dataframe per cell out of your mesh, spots, object or tiff dataset.
#'
#' \code{onePerCell()} returns a dataframe with one table per cell, simplified especially for people who wish
#' to make their own plots (e.g. histograms, etc) based on the summarized cell data. Automatically, it groups
#' the input dataset by unique cell (cell & frame), but it is possible to use this function with other column names.
#'
#' @param inputframe the dataset (in form of a dataframe) to be summarized
#' @param by the grouping value. default = c("cell", "frame")
#' @return the function will return a dataframe (tibble) with the data summarized by the grouping value \code{by}.
#' it will do this by returning the common values per cell (mean intensity, etc), removing values that are useless
#' in one-column context (x/y coordinates of outlines), and calculating the mean, median and standard deviation of
#' numerical values it doesn't recognize. if there are character values with different values per grouped variable,
#' the function will take ONLY THE FIRST VARIABLE.
#'
#' @examples
#' \dontrun{
#' df <- data.frame("cell" = c(rep(1, 5), rep(2, 5)), "frame" = rep(1, 10), "channel" = "GFP", "intensity" = c(1, 3, 4, 5, 6, 1, 2, 7, 3, 10))
#'
#' onePerCell(df)
#' }
onePerCell <- function(inputframe, by = c("cell", "frame")) {

  # remove columns with x/y coordinates
  cols <- colnames(inputframe)
  cols <- cols[!cols %in% c(
    "X",
    "Y",
    "ob_out_x",
    "ob_out_y",
    "num",
    "obpath",
    "X_rot",
    "Y_rot",
    "Xrot_micron",
    "Yrot_micron"
  )]

  inputframe <- inputframe %>% dplyr::select(cols)

  for (n in by) {
    inputframe <- inputframe %>% dplyr::group_by(.data[[n]], add = TRUE)
  }



  # summarize
  inputframe_num <- inputframe %>%
    dplyr::summarize_if(is.double, list(mean = mean, sd = stats::sd)) %>%
    dplyr::ungroup()
  inputframe_num <- inputframe_num[, colSums(inputframe_num != 0) > 0]

  inputframe_char <- inputframe %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::summarize_if(is.character, dplyr::funs(dplyr::first))
  inputframe_out <- dplyr::full_join(inputframe_num, inputframe_char) %>% dplyr::ungroup()

  # check for ambiguous stuff
  chcar <- inputframe_char %>%
    dplyr::ungroup() %>%
    dplyr::select_if(is.character) %>%
    colnames()
  message(paste("Character columns : '", chcar, "' are not used to subset this dataset and instead only the first occuring
                variable of each column per group was kept. To subset by these columns, add them to the input parameter 'by'.",
    collapse = " "
  ))
  return(inputframe_out)
}


#' @export
#'
#'
checkVersionCompatible <- function(oldDataFrame, returnDataFrame = TRUE) {
  if ("Xrotum" %in% colnames(oldDataFrame)) {
    colnames(oldDataFrame)[colnames(oldDataFrame) == "Xrotum"] <- "Xrot_micron"
    colnames(oldDataFrame)[colnames(oldDataFrame) == "Yrotum"] <- "Yrot_micron"
    message("Found old variable names 'Xrotum' and 'Yrotum'.")
    if (returnDataFrame == TRUE) {
      message("Changed old variables for new variables 'Xrot_micron' and 'Yrot_micron'")
      return(oldDataFrame)
    }
  }
  if ("Xrotum" %!in% colnames(oldDataFrame)) {
    message("No compatibility problems found.")
    if (returnDataFrame == TRUE) {
      return(oldDataFrame)
    }
  }
}


## Set the pixel to um conversion
#' @export
addPixels2um <- function(pixelName, pixels2um) {
  if (missing(pixels2um)) {
    pixels2um <- readline("Give the conversion factor from pixel to um:  ")
  }
  if (missing(pixelName)) {
    pixelName <- readline("Give the name of your conversion factor um:  ")
  }

  newp <- c(as.numeric(pixels2um))
  names(newp) <- pixelName
  newML <- append(get(magnificationList, envir = magEnv), newp)
  assign(magnificationList, newML, envir = magEnv)
  print("Currently loaded magnification converters:")
  print(get(magnificationList, envir = magEnv))
}

#' @export
getPixels2um <- function() {
  print("Currently loaded magnification converters:")
  print(get(magnificationList, envir = magEnv))
}



## merge spotfiles with only raw coordinates with mesh file with only raw data. add mesh length/width while on it.
#' @export
spotsInBox <- function(spotdata, meshdata, Xs = "x", Ys = "y", Xm = "X", Ym = "Y", meshInOutput = FALSE) {
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {
      utils::install.packages("shotGroups")
    } else {
      stop("Canceled")
    }
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    inp <- readline("Package 'sp' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {
      utils::install.packages("sp")
    } else {
      stop("Canceled")
    }
  }

  # rewrite colnames if not the same as suggested
  pb <- utils::txtProgressBar(min = 0, max = 100, title = "Total Progress SpotsInBox:")
  if (Xs != "x") {
    colnames(spotdata)[colnames(spotdata) == Xs] <- "x"
  }
  if (Ys != "y") {
    colnames(spotdata)[colnames(spotdata) == Ys] <- "y"
  }
  if (Xm != "X") {
    colnames(meshdata)[colnames(meshdata) == Xm] <- "X"
  }
  if (Ym != "Y") {
    colnames(meshdata)[colnames(meshdata) == Ym] <- "Y"
  }

  if ("max.width" %in% colnames(meshdata) == T) {
    u <- 1
  }
  if ("max.width" %in% colnames(meshdata) == F) {
    u <- 2
  }

  if ("length" %in% colnames(meshdata) == T) {
    a <- 1
  }
  if ("length" %in% colnames(meshdata) == F) {
    a <- 2
  } # if length and max width are already defined, don't touch them.
  utils::setTxtProgressBar(pb, 5)

  # for-loop to re-write
  # cellframe_meshdata <- paste(meshdata$cell, meshdata$frame, sep="_")
  # o <- meshdata %>%
  #  dplyr::group_by(.data$frame, .data$cell) %>%
  #   getSpotsInBox(spotdatap = spotdata[spotdata$frame==meshdata$frame,],
  #           u=u,
  #             a=a,
  #           returnMESH = FALSE)


  if (!requireNamespace("pbapply", quietly = TRUE)) {
    o <- lapply(unique(meshdata$frame), function(x) {
      lapply(unique(meshdata[meshdata$frame == x, ]$cell), function(y) {
        getSpotsInBox(
          meshp = meshdata[meshdata$frame == x & meshdata$cell == y, ],
          spotdatap = spotdata[spotdata$frame == x, ],
          u,
          a,
          returnMESH = meshInOutput
        )
      })
    })
  }
  if (requireNamespace("pbapply", quietly = TRUE)) {
    o <- pbapply::pblapply(unique(meshdata$frame), function(x) {
      lapply(unique(meshdata[meshdata$frame == x, ]$cell), function(y) {
        getSpotsInBox(
          meshp = meshdata[meshdata$frame == x & meshdata$cell == y, ],
          spotdatap = spotdata[spotdata$frame == x, ],
          u,
          a,
          returnMESH = meshInOutput
        )
      })
    })
  }

  utils::setTxtProgressBar(pb, 45)

  outs <- list()
  outs$spots_relative <- do.call("rbind", lapply(o, function(x) do.call("rbind", lapply(x, function(y) y$REP))))
  utils::setTxtProgressBar(pb, 65)
  names(outs) <- c("spots_relative")

  if (meshInOutput == TRUE) {
    outs$mesh <- do.call("rbind", lapply(o, function(x) do.call("rbind", lapply(x, function(y) y$MESH))))
    utils::setTxtProgressBar(pb, 85)
    names(outs) <- c("spots_relative", "mesh")
  }
  utils::setTxtProgressBar(pb, 100)
  return(outs) # return datasets as list of dataframes
}

getSpotsInBox <- function(meshp, spotdatap, u, a, returnMESH = FALSE) {
  box <- suppressWarnings(shotGroups::getMinBBox(data.frame(x = meshp$X, y = meshp$Y))) # bounding box of cell
  lengthwidth <- c(box$width, box$height)

  if (u == 2) {
    meshp$max.width <- min(lengthwidth)
  }
  if (a == 2) {
    meshp$max.length <- max(lengthwidth) # take length/width if not already defined
  }

  pinps <- suppressWarnings(sp::point.in.polygon(spotdatap$x, spotdatap$y, meshp$X, meshp$Y)) # find spot/object coordinates inside cell
  pinps <- data.frame("x" = spotdatap$x, "y" = spotdatap$y, "pip" = pinps)
  if (nrow(pinps) > 0) {
    pinps <- pinps[pinps$pip == 1, ]
  }

  pts <- data.frame(box$pts) # get midpoint of the bounding box + define median lines
  shadowpts <- rbind(pts[2:4, ], pts[1, ])
  distances <- (pts + shadowpts) / 2 # coordinates of median lines

  d1 <- data.frame(x = abs(distances$x - distances$x[1]), y = abs(distances$y - distances$y[1]))
  d1$pointn <- 1:4
  d1 <- d1[-1, ]
  d1$dist <- polar_distance(d1$x, d1$y)


  un <- data.frame(table(round(d1$dist, 5)))

  partdist <- un$Var1[un$Freq == 1]
  partner <- d1$pointn[round(d1$dist, 5) == partdist]
  rest <- d1$pointn[round(d1$dist, 5) != partdist]

  if (round(d1$dist[d1$pointn == partner], 5) == round(min(lengthwidth), 5)) {
    widthline <- distances[c(1, partner), ]
    lengthline <- distances[rest, ]
  }
  if (round(d1$dist[d1$pointn == partner], 5) == round(max(lengthwidth), 5)) {
    widthline <- distances[rest, ]
    lengthline <- distances[c(1, partner), ] # pick which line is length/width
  }

  mp <- c(mean(lengthline$x), mean(lengthline$y)) # midpoint
  X_cor <- meshp$X - mp[1]
  Y_cor <- meshp$Y - mp[2]
  angle <- (-box$angle) * pi / 180 # angle to lay cell flat on x axis

  rotations <- rotate(X_cor, Y_cor, angle)
  meshp$X_rot <- rotations$x
  meshp$Y_rot <- rotations$y

  if (nrow(pinps) > 0) { # rotate spot/object points
    Lc <- pinps$x - mp[1]
    Dc <- pinps$y - mp[2]
    rotations <- rotate(Lc, Dc, angle)
    pinps$l <- -rotations$x
    pinps$d <- -rotations$y
    pinps$max.width <- unique(meshp$max.width)
    if ("max_length" %in% colnames(meshp)) {
      pinps$max.length <- unique(meshp$max_length)
      mesh$max.length <- mesh$max_length
      mesh$max_length <- NULL
    } else {
      pinps$max.length <- unique(meshp$max.length)
    }
    pinps$cell <- unique(meshp$cell)
    pinps$frame <- unique(meshp$frame)
  }
  # if(i==min.i&&n==min.n){ #first data frame in dataset
  #  if(nrow(pinps)>0){
  #    REP <- pinps
  #    }
  #  Mfull <- meshp
  #  }
  if ("trajectory" %in% colnames(spotdatap) == T) {
    pinps <- merge(spotdatap[, c("x", "y", "trajectory", "displacement_sq", "trajectory_length")], pinps)
  }

  if (returnMESH == TRUE) {
    return(list("REP" = pinps, "MESH" = meshp))
  }
  if (returnMESH == FALSE) {
    return(list("REP" = pinps))
  }
}


getpointsaround <- function(datsaround, angle) {
  xlist <- c()
  ylist <- c()

  xlist[1] <- rotate(datsaround$X_corRM, datsaround$Y_corRM, angle)$x
  xlist[2] <- rotate(datsaround$X_corRM, datsaround$Y_corRMa, angle)$x
  xlist[3] <- rotate(datsaround$X_corRMa, datsaround$Y_corRMa, angle)$x
  xlist[4] <- rotate(datsaround$X_corRMa, datsaround$Y_corRM, angle)$x

  ylist[1] <- rotate(datsaround$X_corRM, datsaround$Y_corRM, angle)$y
  ylist[2] <- rotate(datsaround$X_corRM, datsaround$Y_corRMa, angle)$y
  ylist[3] <- rotate(datsaround$X_corRMa, datsaround$Y_corRMa, angle)$y
  ylist[4] <- rotate(datsaround$X_corRMa, datsaround$Y_corRM, angle)$y

  datforpoint <- data.frame(xt = xlist, yt = ylist, pointN = datsaround$pointN)
  return(datforpoint)
}

turnraws <- function(rawdatafile, i, n, mp, angle) {
  rawr <- rawdatafile[rawdatafile$frame == i & rawdatafile$cell == n, ]
  X_corR <- rawr$x - mp[1]
  Y_corR <- rawr$y - mp[2]
  rotations <- rotate(X_corR, Y_corR, angle)
  rawr$X_rot <- rotations$x
  rawr$Y_rot <- rotations$y
  rawr$pointN <- c(1:nrow(rawr))
  datsaround <- data.frame(X_corRM = X_corR - 0.5, X_corRMa = X_corR + 0.5, Y_corRM = Y_corR - 0.5, Y_corRMa = Y_corR + 0.5, pointN = rawr$pointN)
  datsaround <- lapply(1:nrow(datsaround), function(x) getpointsaround(datsaround[x, ], angle))
  datsaround <- do.call(rbind, datsaround)
  rawr <- merge(rawr, datsaround)
  return(rawr)
}

turncell <- function(MESHp, u, rawdatafile, a, n, i, ars) {
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {
      utils::install.packages("shotGroups")
    } else {
      stop("Canceled")
    }
  }
  if (!requireNamespace("sp", quietly = TRUE)) {
    inp <- readline("Package 'sp' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {
      utils::install.packages("sp")
    } else {
      stop("Canceled")
    }
  }
  box <- suppressWarnings(shotGroups::getMinBBox(data.frame(x = MESHp$X, y = MESHp$Y))) # bounding box of cell
  lengthwidth <- c(box$width, box$height)
  if (ars == 2) {
    MESHp$area <- sp::Polygon(cbind(x = MESHp$X, y = MESHp$Y))@area
  }
  if (u == 2) {
    MESHp$max.width <- min(lengthwidth)
  }
  if (a == 1) {
    MESHp$max.length <- max(lengthwidth) # take length/width if not already defined
  }
  pts <- data.frame(box$pts) # get midpoint of the bounding box + define median lines
  shadowpts <- rbind(pts[2:4, ], pts[1, ])
  distances <- (pts + shadowpts) / 2 # coordinates of median lines

  d1 <- data.frame(x = abs(distances$x - distances$x[1]), y = abs(distances$y - distances$y[1]))
  d1$pointn <- 1:4
  d1 <- d1[-1, ]
  d1$dist <- polar_distance(d1$x, d1$y)


  un <- data.frame(table(round(d1$dist, 5)))

  partdist <- un$Var1[un$Freq == 1]
  partner <- d1$pointn[round(d1$dist, 5) == partdist]
  rest <- d1$pointn[round(d1$dist, 5) != partdist]

  if (round(d1$dist[d1$pointn == partner], 5) == round(min(lengthwidth), 5)) {
    widthline <- distances[c(1, partner), ]
    lengthline <- distances[rest, ]
  }
  if (round(d1$dist[d1$pointn == partner], 5) == round(max(lengthwidth), 5)) {
    widthline <- distances[rest, ]
    lengthline <- distances[c(1, partner), ] # pick which line is length/width
  }

  mp <- c(mean(lengthline$x), mean(lengthline$y)) # midpoint
  X_cor <- MESHp$X - mp[1]
  Y_cor <- MESHp$Y - mp[2]
  angle <- (180 - box$angle) * pi / 180 # angle to lay cell flat on x axis

  MESHp$angle <- angle
  MESHp$Xmid <- mp[1]
  MESHp$Ymid <- mp[2]
  rotations <- rotate(X_cor, Y_cor, angle)
  MESHp$X_rot <- rotations$x
  MESHp$Y_rot <- rotations$y
  if (MESHp$Y_rot[[1]] > 0) {
    rotations <- rotate(MESHp$X_rot, MESHp$Y_rot, pi)
    MESHp$Y_rot <- rotations$y
    MESHp$X_rot <- rotations$x
    MESHp$angle <- MESHp$angle - pi
  }


  if (!missing(rawdatafile)) {
    message(paste("Turning raw data for cell", n))

    rawr <- turnraws(rawdatafile, i, n, mp, angle)
    return(list(mesh = MESHp, rawdat = rawr))
  }
  if (missing(rawdatafile)) {
    return(MESHp)
  }
}


meshTurn <- function(MESH, Xm = "X", Ym = "Y", rawdatafile) {
  if (Xm != "X") {
    colnames(MESH)[colnames(MESH) == Xm] <- "X"
  }
  if (Ym != "Y") {
    colnames(MESH)[colnames(MESH) == Ym] <- "Y"
  }

  if ("max.width" %in% colnames(MESH)) {
    u <- 1
  } else {
    u <- 2
  }
  if ("max.length" %in% colnames(MESH)) {
    a <- 2
  } else {
    a <- 1
  }
  if ("area" %in% colnames(MESH)) {
    ars <- 1
  } else {
    ars <- 2
  }
  if ("x0" %in% colnames(MESH)) {
    MESH <- spotrXYMESH(MESH)
  }
  if (!missing(rawdatafile)) {
    Rlist <- list()
  }
  Mlist <- list()

  # if length and max width are already defined, don't touch them.
  for (i in unique(MESH$frame)) { # per frame
    # print(paste("Turning meshes for frame", i))
    M <- MESH[MESH$frame == i, ]
    if (!missing(rawdatafile)) {
      Mlistboth <- lapply(unique(M$cell), function(x) turncell(M[M$cell == x, ], u, rawdatafile, a, x, i, ars = ars))
      MlistF <- lapply(Mlistboth, function(x) x$mesh)
      RlistF <- lapply(Mlistboth, function(x) x$rawdat)
      Rlist[[i]] <- do.call(rbind, RlistF)
    }
    if (missing(rawdatafile)) {
      MlistF <- lapply(unique(M$cell), function(x) turncell(M[M$cell == x, ], u, a = a, n = x, i = i, ars = ars))
    }
    Mlist[[i]] <- do.call(rbind, MlistF)
  }

  Mfull <- do.call(rbind, Mlist)
  if (missing(rawdatafile)) {
    return(Mfull) # return datasets as list of dataframes
  }
  rawdata_turned <- do.call(rbind, Rlist)
  rawdata_turned <- merge(rawdata_turned, unique(Mfull[, c("cell", "frame", "max.length", "max.width", "area")]))

  return(list(mesh = Mfull, rawdata_turned = rawdata_turned))
}

spotrXYMESH <- function(MESH, x_1 = "x1", y_1 = "y1", x_0 = "x0", y_0 = "y0") {
  u <- colnames(MESH)
  MESH <- MESH[!is.na(MESH$cell), ]
  MESH0 <- MESH[, u[u != x_1 & u != y_1]]
  MESH1 <- MESH[, u[u != x_0 & u != y_0]]
  MESH1$xy <- 1
  MESH0$xy <- 0
  colnames(MESH1) <- gsub("1", "", colnames(MESH1))
  colnames(MESH0) <- gsub("0", "", colnames(MESH0))
  ## need to fix lapply function - check dapply? mapply? otherwise first make list of each row.
  MESH1$cf <- paste(MESH1$cell, MESH1$frame, sep = "_")
  MESH1 <- do.call("rbind", lapply(unique(MESH1$cf), function(x) mesh1Fun(MESH1[MESH1$cf == x, ])))
  MESH1$cf <- NULL
  MESH <- merge(MESH0, MESH1, all = T)
  colnames(MESH)[colnames(MESH) == "x"] <- "X"
  colnames(MESH)[colnames(MESH) == "y"] <- "Y"
  MESH <- MESH[order(MESH$frame, MESH$cell, MESH$num), ]
  return(MESH)
}

mesh1Fun <- function(MESH1) {
  MESH1$num <- 1 + max(MESH1$num, na.rm = T) * 2 - MESH1$num
  return(MESH1)
}

mergeframes <- function(REP, MESH, mag = "No_PixelCorrection", cutoff = T, maxfactor = 2, minfactor = 0.5, remOut = T, ouf = F) {

  # REP<- REP[(0<REP$l),]
  if ("rel.l" %in% colnames(REP)) {
    REP <- REP[(REP$rel.l < 1), ]
  }
  if ("max.length" %in% colnames(MESH) != T & "max_length" %in% colnames(MESH) != T) {
    MESH$max.length <- MESH$length
  }
  if ("max_length" %in% colnames(MESH)) {
    MESH$max.length <- MESH$max_length
  }
  if ("area" %in% colnames(MESH)) {
    M <- unique(MESH[, c("cell", "frame", "max.length", "max.width", "area")])
  } else {
    M <- unique(MESH[, c("cell", "frame", "max.length", "max.width")])
  }
  M <- M[order(M$max.length), ]
  M$cellnum <- c(1:nrow(M))

  # merging
  MR <- merge(M, REP, all = T)

  # remove MR's cells which have NA's in the cell area
  MR <- MR[!is.na(MR$max.length), ]
  # remove duplicated rows
  MR <- MR[!duplicated(MR$l) | is.na(MR$l) & !duplicated(MR$d) | is.na(MR$d), ]

  MR <- spotMR(MR)
  MR$totalspot[is.na(MR$totalspot)] <- 0

  # if needed: remove smallest and largest ones (cutoff: smaller than 1/2 mean and larger than 2x mean)
  if (cutoff == T) {
    MR <- MR[MR$max.length < (maxfactor * mean(MR$max.length)), ]
    MR <- MR[MR$max.length > (minfactor * mean(MR$max.length)), ]
  }

  MR <- MR[order(MR$max.length), ]

  # make column with row numbers per cell length.
  MR$num <- c(1:nrow(MR))

  pix2um <- unlist(get(magnificationList, envir = magEnv)[mag])
  MR <- LimDum(MR, pix2um, ouf = ouf)
  return(MR)
}

# add spotnumbers!
spotMR <- function(dat) {
  if ("spot" %in% colnames(dat)) {
    # NA in spots replaced by "0"
    dat$spot[is.na(dat$spot)] <- 0
  } else {
    dat <- dat[order(dat$frame, dat$cell, dat$max.length), ]
    dat$spot <- 1
    for (n in 1:(nrow(dat) - 1)) {
      if (dat$max.length[n + 1] == dat$max.length[n]) {
        dat$spot[n + 1] <- dat$spot[n] + 1
      }
    }
    dat$spot[is.na(dat$spot)] <- 0
  }
  if ("totalspot" %in% colnames(dat)) {
    dat$totalspot[is.na(dat$totalspot)] <- 0
  } else {
    dat$cellframe <- paste(dat$cell, dat$frame, sep = ".")
    spotn <- data.frame(cellframe = unique(dat$cellframe), totalspot = unlist(lapply(unique(dat$cellframe), function(x) max(dat$spot[dat$cellframe == x]))))
    dat <- merge(dat, spotn)
  }
  return(dat)
}

centrefun <- function(dat, xie = "ob_x", yie = "ob_y") {
  if (!requireNamespace("shotGroups", quietly = TRUE)) {
    inp <- readline("Package 'shotGroups' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {
      utils::install.packages("shotGroups")
    } else {
      stop("Canceled")
    }
  }
  dat <- dat[!is.na(dat$ob_x), ]
  dat$centre_x <- NA
  dat$centre_y <- NA
  if (requireNamespace("pbapply", quietly = TRUE)) {
    dat <- do.call("rbind", pbapply::pblapply(unique(dat$obID), function(x) takeObjectCentre(dat[dat$obID == x, ], xie, yie)))
  }
  if (!requireNamespace("pbapply", quietly = TRUE)) {
    dat <- do.call("rbind", lapply(unique(dat$obID), function(x) takeObjectCentre(dat[dat$obID == x, ], xie, yie)))
  }
  return(dat)
}

takeObjectCentre <- function(dat, xie, yie) {
  woei <- as.matrix(dat[c(xie, yie)])
  cen <- shotGroups::getMinCircle(woei)$ctr
  dat$centre_x <- cen[1]
  dat$centre_y <- cen[2]
  return(dat)
}

# add object centre to mesh file and turn accordingly
#' @importFrom dplyr %>%
midobject <- function(MESH, OBJ, p2um) {
  OBJ$angle <- NULL
  OBJ$max.length <- NULL
  MESH <- MESH %>%
    dplyr::left_join(OBJ) %>%
    dplyr::mutate(
      xccor = .data$centre_x - .data$Xmid,
      yccor = .data$centre_y - .data$Ymid,
      xcorcor = .data$ob_x - .data$Xmid,
      ycorcor = .data$ob_y - .data$Ymid,
      Lmid = .data$xccor * cos(.data$angle) - .data$yccor * sin(.data$angle),
      Dum = .data$xccor * sin(.data$angle) + .data$yccor * cos(.data$angle),
      ob_out_x = .data$xcorcor * cos(.data$angle) - .data$ycorcor * sin(.data$angle),
      ob_out_y = .data$xcorcor * sin(.data$angle) + .data$ycorcor * cos(.data$angle)
    ) %>%
    dplyr::select(.data$frame, .data$cell, .data$obpath, .data$obnum, .data$obID, .data$max.length, .data$max.width, .data$Dum, .data$Lmid, .data$ob_out_x, .data$ob_out_y) %>%
    dplyr::distinct() %>%
    dplyr::mutate(num = dplyr::dense_rank(.data$max.length)) %>%
    LimDum(p2um) %>%
    dplyr::mutate(ob_out_x = .data$ob_out_x * p2um, ob_out_y = .data$ob_out_y * p2um)
  return(MESH)

  # MESH$xccor <- MESH$centre_x - MESH$Xmid
  # MESH$yccor <- MESH$centre_y - MESH$Ymid
  # MESH$xcorcor <- MESH$ob_x - MESH$Xmid
  # MESH$ycorcor <- MESH$ob_y - MESH$Ymid

  # MESH$Lmid <- MESH$xccor * cos(MESH$angle) - MESH$yccor * sin(MESH$angle)
  # MESH$Dum <- MESH$xccor * sin(MESH$angle) + MESH$yccor * cos(MESH$angle)
  # MESH$ob_out_x <- MESH$xcorcor * cos(MESH$angle) - MESH$ycorcor * sin(MESH$angle)
  # MESH$ob_out_y <- MESH$xcorcor * sin(MESH$angle) + MESH$ycorcor * cos(MESH$angle)
  # MO <- MESH[,c("frame", "cell", "obpath", "obnum", "obID", "max.length", "max.width", "Dum", "Lmid", "ob_out_x", "ob_out_y")]
  #  MO <- unique(MO)

  # MOnum <- unique(MO[,c("frame", "cell", "max.length", "obnum")])
  # MOnum <- MOnum[order(MOnum$max.length),]
  # MOnum$num <- 1:nrow(MOnum)
  # MO <- merge(MOnum, MO)
  # MO <- LimDum(MO, p2um)
  # MO$ob_out_x <- MO$ob_out_x*p2um
  # MO$ob_out_y <- MO$ob_out_y*p2um
  # return(MO)
}

################################################################################################
# plot preparation
# quartiles, maxima, etc.
LimDum <- function(MR, pix2um, remOut = T, ouf = F) {
  if ("l" %in% colnames(MR) == TRUE & ouf == F) {
    MR$Lmid <- MR$l * pix2um
  }
  if ("l" %in% colnames(MR) == TRUE & ouf == T) {
    MR$Lmid <- (MR$l - (MR$max.length / 2)) * pix2um
  }
  if ("l" %in% colnames(MR) != TRUE) {
    MR$Lmid <- MR$Lmid * pix2um
  }
  MR$pole1 <- -MR$max.length * 0.5 * pix2um
  MR$pole2 <- -MR$pole1
  if ("d" %in% colnames(MR)) {
    MR$Dum <- MR$d * pix2um
  }
  if ("d" %in% colnames(MR) == FALSE) {
    MR$Dum <- MR$Dum * pix2um
  }
  MR$max_um <- MR$max.length * pix2um
  MR$maxwum <- MR$max.width * pix2um
  if (remOut == T) {
    MR <- MR[abs(MR$Lmid) < MR$pole2 | is.na(MR$Lmid), ]
    MR <- MR[abs(MR$Dum) < (MR$max.width / 2) | is.na(MR$Lmid), ]
  }
  return(MR)
}

mL <- list("100x_TIRF" = 0.05204891, "100x_DVMolgen" = 0.0645500, "No_PixelCorrection" = 1, "100x_FRAP" = 0.0499548)
magnificationList <- "magnificationList"
magEnv <- new.env()
assign(magnificationList, mL, envir = magEnv)



#' @export
micron <- function() {
  return("\u00b5m")
}


#' @export
#' @title Orient your cells by the side in which your spots or objects (predominantly) are located
#'
#' \code{orientCells()} takes a \code{spots_relative} or \code{object_relative} dataset and uses the relative localization of each spot or mid-point of each object
#' to flip all related cells & content to one side based on this localization.
#'
#' @param dataTurn has to be a \code{spots_relative} or \code{object_relative} dataset
#' @param cellData a list of corresponding datasets that need to be turned as well. can be 1 dataset (for instance \code{mesh}) or a list of datasets (like \code{list(mesh, object_relative)})
#' @param removeNonOriented default is FALSE - if set to TRUE, in all datasets only the oriented cells will be kept.
#' @return the function will return a list of the datasets, where the relative X/Y coordinates of the cells and internal features are all oriented to one side.
#'
#' @examples
#' \dontrun{
#'
#' output <- orientCells(dataTurn = myData$spots_relative, cellData = list(myData$mesh, myData$object_relative))
#' }
#'
orientCells <- function(dataTurn, cellData, removeNonOriented = FALSE) {
  # check if the dataset is relative
  if (!"Lmid" %in% colnames(dataTurn)) {
    if (!"l" %in% colnames(dataTurn)) {
      stop("No relative localizations found in the 'dataTurn' input. Please use a dataframe with relative spot or object locatizations.")
    }
    if ("l" %>% colnames(dataTurn)) {
      # rename l & d
      dataTurn <- dataTurn %>%
        dplyr::rename(Lmid = .data$l, Dum = .data$d)
    }
  }
  # in case of objects
  if ("obID" %in% colnames(dataTurn)) {
    dataTurn <- dataTurn %>%
      dplyr::group_by(.data$frame, .data$cell) %>%
      dplyr::mutate(
        totalspot = dplyr::n_distinct(.data$obID),
        spot = .data$obID
      ) %>%
      dplyr::ungroup()
  }

  # classify location of spot/object middle with variable 'pol' (handy if you have more than 1 spot)
  dataTurn <- dataTurn %>%
    dplyr::mutate(pol = dplyr::if_else(.data$Lmid < 0, -1, 1, 0))

  if (max(dataTurn$totalspot) > 1) {
    message("More than one spot per cell detected. Function takes the side most spots are on (in case of a tie, cell is classified as 'non-polarized' by setting variable 'pol' to 0).")
    dataSmall <- dataTurn %>%
      dplyr::distinct(.data$frame, .data$cell, .data$spot, .data$pol) %>%
      dplyr::group_by(.data$frame, .data$cell) %>%
      dplyr::mutate(pol = mean(.data$pol)) %>%
      dplyr::mutate(pol = dplyr::if_else(.data$pol == -1 | .data$pol == 1, .data$pol, 0, 0))
    dataTurn <- dataTurn %>%
      dplyr::select(-.data$pol) %>%
      dplyr::left_join(dataSmall)
  }

  if (removeNonOriented == TRUE) {
    dataTurn <- dataTurn[dataTurn$pol != 0, ]
  }

  listPol <- dataTurn %>%
    dplyr::select(.data$frame, .data$cell, .data$pol)

  dataTurn <- turnEach(dataTurn)

  if (is.data.frame(cellData) == FALSE) {
    if (is.list(cellData) == FALSE) {
      stop("Does not recognize the cellData input. Either provide one dataframe or a list() of dataframes.")
    }
    cellData <- lapply(cellData, function(x) turnEach(x, listPol))
  } else {
    cellData <- turnEach(cellData, listPol)
  }

  return(list(dataTurn = dataTurn, cellData = cellData))
}


turnEach <- function(partCD, listPol) {
  if ("spot" %in% colnames(partCD)) {
    dType <- "spot"
    if (!"Lmid" %in% colnames(partCD)) {
      if (!"l" %in% colnames(partCD)) {
        stop("The spot dataframe does not contain relative localizations. Please add a spots_relative dataframe.")
      }
      partCD <- partCD %>%
        dplyr::rename(Lmid = .data$l, Dum = .data$d)
    }
  }
  if ("obID" %in% colnames(partCD)) {
    dType <- "object"
  }
  if ("Xmid" %in% colnames(partCD)) {
    dType <- "mesh"
  }
  if (missing(dType)) {
    stop("Did not recognize the input dataframe in cellData.")
  }

  if (!missing(listPol)) {
    partCD <- partCD %>%
      dplyr::right_join(listPol)
  }


  if (dType == "spot") {
    partCD <- partCD %>%
      dplyr::mutate(
        Lmid = dplyr::if_else(.data$pol != 0, .data$Lmid * .data$pol, .data$Lmid, .data$Lmid),
      )
  }
  if (dType == "object") {
    partCD <- partCD %>%
      dplyr::mutate(
        Lmid = dplyr::if_else(.data$pol != 0, .data$Lmid * .data$pol, .data$Lmid, .data$Lmid),
        ob_out_x = dplyr::if_else(.data$pol != 0, .data$ob_out_x * .data$pol, .data$ob_out_x, .data$ob_out_x)
      )
  }
  if (dType == "mesh") {
    partCD <- partCD %>%
      dplyr::mutate(
        X_rot = dplyr::if_else(.data$pol != 0, .data$X_rot * .data$pol, .data$X_rot, .data$X_rot),
        Xrot_micron = dplyr::if_else(.data$pol != 0, .data$Xrot_micron * .data$pol, .data$Xrot_micron, .data$Xrot_micron)
      )
    if ("xt" %in% colnames(partCD)) {
      partCD <- partCD %>%
        dplyr::mutate(
          xt = dplyr::if_else(.data$pol != 0, .data$xt * .data$pol, .data$xt, .data$xt)
        )
    }
  }

  return(partCD)
}






polar_distance <- function(x, y) {
  return(sqrt(x^2 + y^2))
}

polar_angle <- function(x, y) {
  return(atan(x / y))
}

rotate <- function(x, y, theta, type = "rad") {
  theta <- switch(type,
    "rad" = theta,
    "deg" = (pi / 180) * theta,
    "hour" = 0.2618 * theta
  )
  rotated_x <- x * cos(theta) - y * sin(theta)
  rotated_y <- x * sin(theta) + y * cos(theta)
  return(data.frame(x = rotated_x, y = rotated_y))
}

st_rotate_around <- function(geometry, theta, around = sf::st_point(c(0, 0))) {
  rotation_matrix <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  rotated_geometry <- ((geometry - around) * rotation_matrix) + around
  # if (sf::st_geometry_type(rotated_geometry) == "POLYGON") {
  #   colnames(rotated_geometry[[1]]) <- c("x", "y")
  # }
  return(rotated_geometry)
}



get_minimum_bounding_box_centroid_and_angle <- function(x, y) {
  bounding_box <- shotGroups::getMinBBox(data.frame(point.x = x, point.y = y))
  points <- as.matrix(bounding_box$pts)
  points <- rbind(points, points[1, ])
  centroid <- sf::st_polygon(list(points)) |> sf::st_centroid()
  return(list(centroid = centroid, angle = bounding_box$angle))
}

st_minimum_bounding_box <- function(polygon) {
  polygon_df <- polygon |> 
    unlist() |>
    matrix(ncol = 2) |>
    as_tibble() |>
    rename(c("point.x" = 1, "point.y" = 2))
  bounding_box_data <- shotGroups::getMinBBox(polygon_df)
  bonding_box_polygon <- st_polygon_autoclose(bounding_box_data$pts[,"x"], bounding_box_data$pts[,"y"])
  return(list(polygon = list(bonding_box_polygon), angle = bounding_box_data$angle))
}







angle_from_horizontal <- function(..., type = "rad") {
  result <- switch(type,
    "rad" = pi,
    "deg" = 180,
    "hour" = 12
  )
  angles <- c(...)
  for (theta in angles) {
    result <- result - theta
  }
  return(result)
}

angle_from_vertical <- function(..., type = "rad") {
  result <- switch(type,
    "rad" = pi / 2,
    "deg" = 90,
    "hour" = 6
  )
  angles <- c(...)
  for (theta in angles) {
    result <- result - theta
  }
  return(result)
}

st_polygon_autoclose <- function(x, y) {
  if (first(x) != last(x) || first(y) != last(y)) {
    x <- append(x, first(x))
    y <- append(y, first(y))
  }
  return(sf::st_polygon(list(cbind(x = x, y = y))))
}

is_polygon_rotation_y_positive <- function(polygon, theta, around = sf::st_point(c(0, 0))) {
  rotated_polygon <- st_rotate_around(polygon, theta, around)
  return(rotated_polygon[[1]][1, "y"] > 0)
}

st_triangulate_circle <- function(center, radius, angle, bins, angle_unit = "rad") {
  if(!missing(angle)) {
    angle <- switch(angle_unit,
    "rad" = angle,
    "deg" = π/180 * angle)
  }

  angle_ <- switch(rlang::check_exclusive(angle, bins),
    angle = angle,
    bins = (2 * pi) / bins)

  bins_ <- switch(rlang::check_exclusive(angle, bins),
    angle = 2 * pi / angle,
    bins = bins)

  top_left_point <- st_rotate_around(center + sf::st_point(c(0, radius)), -angle_ / 2, around = center)

  splits <- purrr::map(1:bins_, \(x) st_rotate_around(top_left_point, x * angle_, around = center))
  angles <- c(1:bins_) * angle_

  triangles <- tibble(point = splits, next_point = lead(splits, default = splits[1]), angle = angles) |>
    rowwise() |>
    mutate(triangle = list(sf::st_polygon(list(matrix(rbind(unlist(c(center, point, next_point, center))), ncol = 2))))) |>
    select(angle, triangle)

  return(st_sf(triangles))
}

get_points_inside_polygons <- function(raster, polygons_sf, join = st_within, left = FALSE) {
  names(raster) <- "value"
  points <- st_as_sf(as.points(raster))
  points_inside <- st_join(points, polygons_sf, join = join, left = left) |> as_tibble() |> rename(point = geometry)
  return(points_inside)
}

correct_position <- function(points_inside, point, angle, center) {
  points_inside |>
    rowwise() |>
    mutate(
      corrected_pos = list(st_rotate_around(st_point({{point}}[[1]]), -{{angle}}, center)),
      x = round(pos_corr[[1]][1]),
      y = round(pos_corr[[1]][2])) |>
    select(-corrected_pos) -> corrected
  return(corrected)
}

st_ring <- function(inner_radius, outer_radius, center = st_point(c(0, 0))) {
  outer <- st_buffer(center, dist = outer_radius, endCapStyle="ROUND")
  inner <- st_buffer(center, dist = inner_radius, endCapStyle="ROUND")
  ring <- st_difference(outer, inner)
  return(ring)
}

st_concentrics_rings <- function(max_radius, n, center =  st_point(c(0, 0))) {
  size <- max_radius/n
  steps <- seq(from = 0, to = max_radius, by = size)
  rings <- tibble(
    inner = steps[1:length(steps) - 1],
    outer = steps[2:length(steps)]) 
  rings |>
    rowwise() |>
    mutate(ring = list(st_ring(inner, outer, center))) |> st_sf() -> rings
  return(rings)
}

st_box <- function(left, right, top, bot, center, width, height, angle = 0) {
  if (!missing(center)) {
    if (!missing(width) & !missing(height)) {
      left <- center[1] - width/2
      right <- center[1] +  width/2
      bot <- center[2] - height/2
      top <- center[2] + height/2
    }
    else if (!missing(left) & !missing(right)) {
      right <- center[1] + (center-left)/2
      top <- center[2] + (center-bot)/2
    }
    else {
      stop("you need to provide 'width' and 'height' OR 'left' and 'bot' arguments when using 'center' argument")
    }
  }
  else {
    
    if (!xor(missing(top) | missing(right), missing(width) | missing(height))) {
      stop("you need to provide 'width' and 'height' OR 'top' arguments when using 'left' and bot argument")
      }
    if (missing(top) | missing(right)) {
      right <- left +  width/2
      top <- bot +  height/2
    }
  }

  points <- matrix(c(right, top, left, top, left, bot, right, bot, right, top), ncol = 2, byrow = TRUE)
  ox <- (right + left)/2
  oy <- (top + bot)/2
  points <- matrix(c(
    cos(-angle) * (points[,1] - ox)  - sin(-angle) * (points[,2] - oy) + ox,
    sin(-angle) * (points[,1] - ox)  + cos(-angle) * (points[,2] - oy) + oy), ncol = 2)
  box <- st_polygon(list(points))
  return(box)
}


# triangles <- st_triangulate_circle(center, 100, bins = 360)
# rings <- st_concentrics_rings(100, 20, center = center)
# triangles_inside <- get_points_inside_polygons(raster, triangles)
# rings_inside <- get_points_inside_polygons(raster, rings)

# merged <- inner_join(triangles_inside, rings_inside, by = (c("point" = "point", "value" = "value")))
# merged |>
#   group_by(angle, inner) |>
#   summarise(value = mean(value)) |>
#   ggplot() +
#     geom_tile(
#       aes(
#         x = angle,
#         y = inner,
#         fill = value)) +
#     scale_fill_viridis_c() +
#     coord_polar()


center_and_orient_mesh <- function(mesh, angle = 0, around, offset) {
  if ("sfc" %in% class(mesh)) {
    is_sfc <- TRUE
    mesh <- mesh[[1]]}
  if (missing(around)) {around <-  sf::st_centroid(mesh)}
  if (missing(offset)) {offset <-  sf::st_centroid(mesh)}
  if ("sfc" %in% class(around)) {around <- around[[1]]}
  if ("sfc" %in% class(offset)) {offset <- offset[[1]]}
  out <- st_rotate_around(mesh, angle, around) - offset
  if (is_sfc) {out <- st_sfc(out)}
  return(out)
}


morpho_out |>
  rowwise() |> 
  mutate(bb =list(st_minimum_bounding_box(mesh_cont))) |> 
  unnest_wider(bb) |> 
  unnest(polygon) |>
  mutate(polygon = st_sfc(polygon)) |>
  rowwise() |> 
  mutate(
    angle = angle * (2*pi/360), 
    mesh_cent = center_and_orient_mesh(mesh_cont, angle)) |>
  st_as_sf() -> morpho



mj_out  |>
  rowwise() |> 
  mutate(bb =list(st_minimum_bounding_box(mesh))) |> 
  unnest_wider(bb) |> 
  unnest(polygon) |>
  rowwise() |> 
  mutate(
    angle = angle * (2*pi/360), 
    mesh_cent = st_sfc(center_and_orient_mesh(mesh, angle,offset = st_centroid(polygon), around = st_centroid(polygon)))) |>
  st_as_sf() -> mj


StatStackMesh <- ggproto("StatStackMesh", Stat,
  compute_group = function(data, scales) {
  print(tibble(data))
  data |>
    rowwise() |> 
    mutate(bb =list(st_minimum_bounding_box(geometry))) |> 
    unnest_wider(bb) |> 
    unnest(polygon) |>
    rowwise() |> 
    mutate(
    angle = angle * (2*pi/360), 
    geometry = st_sfc(center_and_orient_mesh(geometry, angle, offset = st_centroid(geometry), around = st_centroid(geometry)))) |>
    select(-polygon) -> data
    print(tibble(data))
    data
  },
  required_aes = c("geometry")
)

geom_stackmesh <- function(mapping = NULL, data = NULL, stat = "sf",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  c(
  layer_sf(
    stat = StatStackMesh, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  ),
  coord_sf(default = TRUE)
  )
}

StatSfStack <- ggproto("StatSfStack", Stat,
  compute_layer = function(self, data, params, layout) {
    # add coord to the params, so it can be forwarded to compute_group()
    params$coord <- layout$coord
    ggproto_parent(Stat, self)$compute_layer(data, params, layout)
  },

  compute_group = function(data, scales, coord) {
    data |> 
      rowwise() |> 
      mutate(bb = list(st_minimum_bounding_box(geometry))) |> 
      unnest_wider(bb) |> 
      select(-polygon) |>
      rowwise() |> 
      mutate(
        angle = angle * (2*pi/360), 
        geometry = st_sfc(center_and_orient_mesh(
          geometry,
          angle, 
          offset = st_centroid(geometry),
          around = st_centroid(geometry)))) -> data
    geometry_data <- data[[ geom_column(data) ]]

    geometry_crs <- sf::st_crs(geometry_data)

    bbox <- sf::st_bbox(geometry_data)

    if (inherits(coord, "CoordSf")) {
      # if the coord derives from CoordSf, then it
      # needs to know about bounding boxes of geometry data
      coord$record_bbox(
        xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
        ymin = bbox[["ymin"]], ymax = bbox[["ymax"]]
      )

      # to represent the location of the geometry in default coordinates,
      # we take the mid-point along each side of the bounding box and
      # backtransform
      bbox_trans <- sf_transform_xy(
        list(
          x = c(rep(0.5*(bbox[["xmin"]] + bbox[["xmax"]]), 2), bbox[["xmin"]], bbox[["xmax"]]),
          y = c(bbox[["ymin"]], bbox[["ymax"]], rep(0.5*(bbox[["ymin"]] + bbox[["ymax"]]), 2))
        ),
        coord$get_default_crs(),
        geometry_crs
      )

      # record as xmin, xmax, ymin, ymax so regular scales
      # have some indication of where shapes lie
      data$xmin <- min(bbox_trans$x)
      data$xmax <- max(bbox_trans$x)
      data$ymin <- min(bbox_trans$y)
      data$ymax <- max(bbox_trans$y)
    } else {
      # for all other coords, we record the full extent of the
      # geometry object
      data$xmin <- bbox[["xmin"]]
      data$xmax <- bbox[["xmax"]]
      data$ymin <- bbox[["ymin"]]
      data$ymax <- bbox[["ymax"]]
    }

    data
  }, 

  required_aes = c("geometry")
)

stat_sf_stack <- function(mapping = NULL, data = NULL, geom = "rect",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
  layer_sf(
    stat = StatSfStack,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


geom_sf_stack <- function(mapping = aes(), data = NULL, stat = "sf_stack",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
  c(
    layer_sf(
      geom = GeomSf,
      data = data,
      mapping = mapping,
      stat = stat,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        ...
      )
    ),
    coord_sf(default = TRUE)
  )
}

geom_column <- function(data) {
  w <- which(vapply(data, inherits, TRUE, what = "sfc"))
  if (length(w) == 0) {
    "geometry" # avoids breaks when objects without geometry list-column are examined
  } else {
    # this may not be best in case more than one geometry list-column is present:
    if (length(w) > 1)
      warn("more than one geometry column present: taking the first")
    w[[1]]
  }
}

morpho_out |>
  rowwise() |> 
  mutate(bb = list(st_minimum_bounding_box(mesh_cont))) |> 
  unnest_wider(bb) |> 
  rowwise() |> 
  mutate(
  angle = angle * (2*pi/360), 
  centroid =  st_centroid(mesh_cont),
  mesh_cent = st_sfc(center_and_orient_mesh(mesh_cont, angle, offset = centroid, around = centroid))) -> morpho_cent


st_join(st_sf(morpho_cent,  sf_column_name = "mesh_cont"), tmp_id, left = FALSE) |> as_tibble() -> inside
inside |> inner_join(tmp_id, by = c("id"= "id")) -> full_data



full_data |> 
    group_by(cellID) |>
    summarize(
      pts = st_sfc(st_multipoint(matrix(unlist(geometry), ncol = 2, byrow = TRUE))),
      order = list(id),
      value = list(lyr.1.y)) |> inner_join(morpho_cent, by = c("cellID" = "cellID")) -> out
    
out |> 
  select(cellID, pts, order, value, mesh_cent, angle, centroid) |>
  rowwise() |>
  mutate(
    pts_cent = st_rotate_around(pts, angle, centroid) - centroid
  ) -> out2
  
out2 |> mutate(pts_cent = cbind(pts_cent))|>
      unnest_longer(c(pts_cent, value)) |>
      rowwise() |>
      mutate(pts_cent = list(st_point(c(pts_cent))),
             box = st_sfc(st_box(center = pts_cent, width = 1, height = 1, angle = angle))) -> out3
  
out3 |> ggplot() +
  geom_sf(aes(geometry = box, fill = value ))

    

full_data |> 
    group_by(cellID) |> 
    mutate(
      center = sf::st_sfc(st_rotate_around(geometry, angle, centroid) - centroid)) -> data_cent

data_cent |> mutate(box = sf::st_sfc(st_box(center = center, width = 1, height = 1))) -> data_cent

data_cent |> mutate(box = sf::st_sfc(st_rotate_around(box, angle, st_centroid(box)))) -> data_cent



data_cent |>
  ggplot() +
    geom_tile(
      aes(
        x = x,
        y = y,
        width = 1,
        height = 1,
        fill = lyr.1.y),
      alpha = 0.3) +
    scale_fill_viridis_c() 

