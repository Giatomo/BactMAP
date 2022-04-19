##8-11-2017

##Renske van Raaphorst

##Extracting Morphometrics Data


##########################################################


##Dependencies:
    #library(R.matlab)


extr_Morphometrics_cellList <- function(morphpath, areacutoff = 24){
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    inp <- readline("Package 'R.matlab' needed for this function to work. Press 'y' to install, or any other key to cancel.")
    if (inp %in% c("y", "Y")) {utils::install.packages("R.matlab")}else{stop("Canceled")}
  }
  morphdata <- R.matlab::readMat(morphpath)
  framenum <- 2 * (1:(length(morphdata$frame)/2))
  for (u in framenum){
    cellListf <- data.frame(t(data.frame(morphdata$frame[u])))
    if (nrow(cellListf) != 0) {
      cellListf$frame <- u/2
      cellList <- ifelse(u == 2, cellListf, rbind(cellList, cellListf))
  }
  cellList <- cellList[cellList$area > areacutoff,]
  return(cellList)
  }
}


spotrExtractMorphMESH <- function(cellList){

  if ("length"%in%colnames(cellList)!=T) {
      MESH <- do.call('rbind', lapply(c(1:nrow(cellList)), function(x) extrMorphCell(cellList[x,], l=FALSE)))
      M1 <- MESH[MESH$numpoint == MESH$pole1,]
      M2 <- MESH[MESH$numpoint == MESH$pole2,]
      M1 <- M1[order(M1$cell ,M1$frame),]
      M2 <- M2[order(M2$cell, M2$frame),]
      M1$max.length <- polar_distance((M1$X-M2$X),  (M1$Y-M2$Y))
      MESH <- merge(MESH, M1[, c("cell", "frame", "max.length")])
  } else {
    MESH <- do.call('rbind', lapply(c(1:nrow(cellList)), function(x) extrMorphCell(cellList[x,], l=TRUE)))
    }
  return(MESH)

}

extrMorphCell <- function(cell, l = TRUE){
  print(c(cell$cellID, cell$frame))
  meshcell <- data.frame(cell$Xcont, cell$Ycont, cell$cellID, cell$area, cell$pole1, cell$pole2, cell$frame)

  #meshcell <- data.frame(cell$Xcont[1], cell$Ycont[1], cell$cellID[1], cell$area[1],
                     #    cell$pole1[1], cell$pole2[1], cell$frame[1])
  colnames(meshcell) <- c("X", "Y", "cell", "area", "pole1", "pole2", "frame")
  if (isTRUE(l)) {
    meshcell$max.length <- cell$length
    meshcell$max.width <- cell$width
  }
  meshcell$numpoint <- 1:nrow(meshcell)
  return(meshcell)
}



#' @export
extr_Morphometrics <- function(morphpath, mag, turncells = TRUE, cellList=FALSE){
  if (!missing(mag)) {
    magnification <- unlist(get(magnificationList, envir=magEnv)[mag])
  }
  if (!is.numeric(mag)) {
    stop("Magnification conversion factor not recognized. Please use addPixels2um('pixelName', pixelsize) to add your conversion factor")
  }
  C <- extr_Morphometrics_cellList(morphpath)
  M <- spotrExtractMorphMESH(C)
  listM <- list()
  if (isTRUE(cellList)) {
    listM$cellList <- C
  }
  if (isTRUE(turncells)) {
    listM$mesh <- meshTurn(M)
    listM$mesh$Xrot_micron <- listM$mesh$X_rot * magnification
    listM$mesh$Yrot_micron <- listM$mesh$Y_rot * magnification
    listM$mesh$maxwum <- listM$mesh$max.width * magnification
  }
  else {
    listM$mesh <- M
  }
  listM$mesh$max_um <- listM$mesh$max.length * magnification
  listM$mesh$area_um <- listM$mesh$area * magnification^2
  listM$pixel2um <- magnification

  listM$mesh$num <- listM$mesh$numpoint
  listM$mesh$numpoint <- NULL
  return(listM)
}




max(unlist(purrr::map(x, \(x) length(x))))

get_max_cell_length_col <- function(col) {
  return(max(unlist(purrr::map(col, \(cell) length(cell)))))
}



purrr::imap(morph_mat_list$frame[2,1,], \(x, i) as_tibble(t(x[,1,])) %>% 
  mutate(frame = i)) %>%
  bind_rows() %>% 
  mutate(across(
    where(\(col) get_max_cell_length_col(col) == 1),
    \(x) unlist(x))) %>% 
  rowwise() %>%
  mutate(across(
    where(\(col) get_max_cell_length_col(col) > 1),
    \(x) list(as.numeric(x)))) -> df


df %>% mutate(mesh = list(tibble(Xcont, Ycont, cellID, area, pole1, pole2, frame))) %>% .$mesh %>% bind_rows() -> meshes

meshes %>% 
  group_by(frame, cellID) %>% 
  mutate(
    numpoint = row_number(),
    x0 = if_else(
      numpoint == pole1,
      Xcont,
      NA_real_),
    x1 = if_else(
      numpoint == pole2,
      Xcont,
      NA_real_),
    y0 = if_else(
      numpoint == pole1,
      Ycont,
      NA_real_),
    y1 = if_else(
      numpoint == pole2,
      Ycont,
      NA_real_)) %>%
    fill(x0, x1, y0, y1, .direction = ("downup")) %>%
    mutate(max.length = polar_distance(x0 - x1, y0-y1)) -> temp

