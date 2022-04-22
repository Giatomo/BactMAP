# plot tracks
# 12.18.18
# R van Raaphorst
# bactMAP

#' @export
plotTracks <- function(meshdata,
                       spotdata,
                       objectdata,
                       tracks = FALSE,
                       ignore_singles = FALSE,
                       movie = FALSE,
                       timepalette_lines = "viridis",
                       timepalette_fill = "magma",
                       dimension = c("length", "width"),
                       turn_cells = TRUE,
                       mag,
                       cell = "all",
                       transparency = 0.2) {

  have_lenght <- "length" %in% dimension
  have_width <- "width" %in% dimension
  have_channel <- "channel" %in% colnames(spotdata)
  have_condition <- "condition" %in% colnames(spotdata)
  have_pole <- "pole1" %in% colnames(spotdata)
  color_mode <- if_else(have_channel, "channel", "frame")


  outplot <- ggplot2::ggplot() +
    ggplot2::theme_minimal()


  if (cell[1] != "all") {
    if (!is.numeric(cell)) {
      stop("Cell number not recognized. Give cell number or keep default (cell='all') to return all cells")
    }
    if (!missing(meshdata)) {
      meshdata <- meshdata[meshdata$cell %in% cell, ]
    }
    if (!missing(objectdata)) {
      objectdata <- objectdata[objectdata$cell %in% cell, ]
    }
    if (!missing(spotdata)) {
      spotdata <- spotdata[spotdata$cell %in% cell, ]
    }
  }


  

  if (!missing(meshdata)) {
    if (have_lenght & have_width) {
      if (turn_cells) {
        # Assume name shoud be "Xrot_micron" & "Yrot_micron" and df already cleaned
        outplot <- outplot +
          ggplot2::geom_polygon(
            data = meshdata,
            ggplot2::aes_string(
              x = "Xrot_micron",
              y = "Yrot_micron",
              color = "frame",
              group = "frame"),
            fill = "black", alpha = transparency / 5
          )
      }
      else {
        outplot <- outplot +
          ggplot2::geom_polygon(
            data = meshdata,
            ggplot2::aes_string(x = "X", y = "Y", color = "frame", group = "frame"),
            fill = "black", alpha = transparency / 5
          )
      }
    }
  }

  if (!missing(objectdata)) {
    if (turn_cells | (have_lenght & !have_width) | (!have_lenght & have_width)) {
      if ("ob_out_x" %!in% colnames(objectdata)) {
        if ("ob_x" %!in% colnames(objectdata)) {
          stop("Cannot find object cobjectdatardinate data (ob_out_x/ob_out_y for turned cells, ob_x/ob_y for raw cobjectdatardinates.")
        }
        if (missing(meshdata)) {
          stop("Cannot find turned object data or mesh information to connect the objects to the cell localizations.
                     Replace object dataframe for a dataframe (object_relative) containing this information or add mesh data to convert dataframe.")
        }
        if (missing(mag)) {
          stop("Need pixel to micron conversion factor 'mag' to correctly connect mesh to object data.")
        }

        objectdata <- suppressWarnings(centrefun(objectdata))
        objectdata <- suppressWarnings(midobject(meshdata, objectdata, get(magnificationList, envir = magEnv)[mag]))
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep = "_")

      if (have_lenght & have_width) {
        outplot <- outplot +
          ggplot2::geom_polygon(
            data = objectdata,
            ggplot2::aes_string(x = "ob_out_x", y = "ob_out_y", group = "frameOB", fill = "frame"),
            color = "black", alpha = transparency
          )
      }
      if (have_lenght & !have_width) {
        Om <- stats::aggregate(objectdata[, c("ob_out_x", "pole1", "pole2")], by = list("cell" = objectdata$cell, "frame" = objectdata$frame, "obID" = objectdata$obID), FUN = max)
        Omin <- stats::aggregate(objectdata[, c("ob_out_x", "pole1", "pole2")], by = list("cell" = objectdata$cell, "frame" = objectdata$frame, "obID" = objectdata$obID), FUN = min)
        Om$obxmax <- Om$ob_out_x
        Om$ob_out_x <- NULL
        Omin$obxmin <- Omin$ob_out_x
        Omin$ob_out_x <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1
        outplot <- outplot +
          ggplot2::geom_rect(
            data = objectdata,
            ggplot2::aes_string(xmin = "obxmin", xmax = "obxmax", ymin = "framemin", ymax = "framemax", group = "frameOB", fill = "frame")
          )
        if (missing(spotdata)) {
          objectdata <- objectdata[order(objectdata$cell, objectdata$frame), ]
          outplot <- outplot +
            ggplot2::geom_path(data = objectdata, ggplot2::aes_string(x = "pole1", y = "frame")) +
            ggplot2::geom_path(data = objectdata, ggplot2::aes_string(x = "pole2", y = "frame"))
          outplot <- outplot + ggplot2::xlab("location on cell length axis (\u00b5m)") + ggplot2::ylab("Time (Frames)")
        }
      }
      if (!have_lenght & have_width) {
        Om <- stats::aggregate(objectdata[, c("ob_out_y", "maxwum")], by = list("cell" = objectdata$cell, "frame" = objectdata$frame, "obID" = objectdata$obID), FUN = max)
        Omin <- stats::aggregate(objectdata[, c("ob_out_y", "maxwum")], by = list("cell" = objectdata$cell, "frame" = objectdata$frame, "obID" = objectdata$obID), FUN = min)
        Om$obymax <- Om$ob_out_y
        Om$ob_out_y <- NULL
        Omin$obymin <- Omin$ob_out_y
        Omin$ob_out_y <- NULL
        objectdata <- merge(Om, Omin)

        objectdata$framemin <- objectdata$frame - 0.1
        objectdata$framemax <- objectdata$frame + 0.1

        outplot <- outplot +
          ggplot2::geom_rect(
            data = objectdata,
            ggplot2::aes_string(xmin = "obymin", xmax = "obymax", ymin = "framemin", ymax = "framemax", group = "frameOB", fill = "frame")
          )

        if (missing(spotdata)) {
          objectdata$pole1 <- objectdata$maxwum / 2
          objectdata$pole2 <- objectdata$pole1 * -1
          objectdata <- objectdata[order(objectdata$cell, objectdata$frame), ]

          outplot <- outplot +
            ggplot2::geom_path(data = objectdata, ggplot2::aes_string(x = "pole1", y = "frame")) +
            ggplot2::geom_path(data = objectdata, ggplot2::aes_string(x = "pole2", y = "frame")) +
            ggplot2::xlab("location on cell length axis (\u00b5m)") +
            ggplot2::ylab("Time (Frames)")
        }
      }
    }
    if (!turn_cells & have_lenght & have_width) {
      if ("ob_x" %in% colnames(objectdata) != TRUE) {
        stop("Cannot find object coordinates 'ob_x'/'ob_y'.")
      }
      objectdata$frameOB <- paste(objectdata$frame, objectdata$obID, sep = "_")
      outplot <- outplot +
        ggplot2::geom_polygon(data = objectdata, ggplot2::aes_string(x = "ob_x", y = "ob_y", group = "frameOB", fill = "frame"), color = "black", alpha = transparency)
    }

    outplot <- outplot + ggplot2::scale_fill_viridis_c(option = timepalette_fill)
  }

  if (!missing(spotdata)) {

    if (tracks) {
      if (ignore_singles) {
        spotdata <- spotdata[spotdata$trajectory != -1, ]
      }
    }
    if (turn_cells | (have_lenght & !have_width) | (!have_lenght & have_width)) {
      if ("Lmid" %in% colnames(spotdata) != T) {

        if (missing(mag)) {
          stop("Please specify the pixel-micron conversion value 'mag'.")
        }

        if ("l" %in% colnames(spotdata) != T) {
          if (missing(meshdata) != T) {
            A <- readline("Data doesn't include relative spot positions. Press 'Y'+enter to start spotsInBox() to relate spot positions to cells or any other key to stop the function.")
            if (A %!in% c("Y", "y")) {
              stop("Function stopped.")
            }
            spotdata <- spotsInBox(spotdata, meshdata)$spots_relative
          }
        }
        spotdata$Lmid <- spotdata$l * unlist(get(magnificationList, envir = magEnv)[mag])
        spotdata$Dum <- spotdata$d * unlist(get(magnificationList, envir = magEnv)[mag])
      }

      if (have_lenght & have_width) {
        outplot <- outplot + ggplot2::geom_point(data = spotdata, ggplot2::aes_string(x = "Lmid", y = "Dum", color = "frame"), size = 2)
        if (tracks == TRUE) {
          outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory != -1, ], ggplot2::aes_string(x = "Lmid", y = "Dum", group = "trajectory", color = "frame"), size = 1)
        }
      }
      if (have_lenght & !have_width) {
          outplot <- outplot + ggplot2::geom_point(data = spotdata, ggplot2::aes_string(x = "Lmid", y = "frame", color = color_mode), size = 2)

        if (tracks) {
          outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory != -1, ], ggplot2::aes_string(x = "Lmid", y = "frame", group = "trajectory", color = color_mode), size = 1)
        }

        if (!have_pole) {
          if (!"maxum" %in% colnames(spotdata)) {
            spotdata$maxum <- spotdata$max.length * unlist(get(magnificationList, envir = magEnv)[mag])
          }
          spotdata$pole1 <- 0.5 * spotdata$maxum
          spotdata$pole2 <- -spotdata$pole1
          spotdata <- spotdata[order(spotdata$cell, spotdata$frame), ]
          outplot <- outplot + ggplot2::geom_path(data = spotdata, ggplot2::aes_string(x = "pole1", y = "frame")) + ggplot2::geom_path(data = spotdata, ggplot2::aes_string(x = "pole2", y = "frame"))
        }
      }

      if (!have_lenght & have_width) {

        outplot <- outplot + ggplot2::geom_point(data = spotdata, ggplot2::aes_string(x = "Dum", y = "frame", color = color_mode), size = 2)

        if (tracks) {
          outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory != -1, ], ggplot2::aes_string(x = "Dum", y = "frame", group = "trajectory", color = color_mode), size = 1)
        }
        if (!have_pole) {
          if (!"maxwum" %in% colnames(spotdata)) {
            spotdata$maxwum <- spotdata$max.width * unlist(get(magnificationList, envir = magEnv)[mag])
          }
          spotdata$pole1 <- 0.5 * spotdata$maxwum
          spotdata$pole2 <- -spotdata$pole1
          spotdata <- spotdata[order(spotdata$cell, spotdata$frame), ]
          outplot <- outplot + ggplot2::geom_path(data = spotdata, ggplot2::aes_string(x = "pole1", y = "frame")) + ggplot2::geom_path(data = spotdata, ggplot2::aes_string(x = "pole2", y = "frame"))
        }
      }
    }
    if (!turn_cells & have_lenght & have_width) {
      outplot <- outplot + ggplot2::geom_point(data = spotdata, ggplot2::aes_string(x = "x", y = "y", color = "frame"), size = 2)
      if (tracks == TRUE) {
        outplot <- outplot + ggplot2::geom_path(data = spotdata[spotdata$trajectory != -1, ], ggplot2::aes_string(x = "x", y = "y", group = "trajectory", color = "frame"), size = 1)
      }
    }
  }

  if (!missing(meshdata) | !missing(spotdata)) {
    if ((have_lenght & !have_width) | (!have_lenght & have_width)) {

      outplot <- if_else(
          have_channel,
          outplot + ggplot2::scale_color_manual(values = colorsCUD()),
          outplot + ggplot2::scale_color_viridis_c(option = timepalette_lines)
      )

    }
    if (have_lenght & have_width) {
      outplot <- outplot + ggplot2::scale_color_viridis_c(option = timepalette_lines)
    }
  }

  if (turn_cells & have_lenght & have_width) {

    outplot <- case_when(
      have_channel  & have_condition  ~ outplot + ggplot2::facet_grid(channel ~ condition ~ cell),
      have_channel  & !have_condition ~ outplot + ggplot2::facet_grid(channel ~ cell),
      !have_channel & have_condition  ~ outplot + ggplot2::facet_grid(condition ~ cell),
      !have_channel & !have_condition ~ outplot + ggplot2::facet_grid(~cell)
    ) + 
    ggplot2::xlab(paste("x (", micron(), ")")) + ggplot2::ylab(paste("y (", micron(), ")")) +
    ggplot2::coord_fixed()
  }

  if (!turn_cells | (have_lenght & !have_width) | (!have_lenght & have_width)) {

    outplot <- if_else(
      have_condition,
      outplot + ggplot2::facet_grid(condition ~ cell, scales = "free"),
      outplot + ggplot2::facet_grid(~cell, scales = "free")
    )

    outplot <- case_when(
      have_lenght  & have_width  ~  outplot + ggplot2::xlab("x (pixels)") + ggplot2::ylab("y (pixels)"),
      have_lenght  & !have_width ~ outplot + ggplot2::xlab(paste("distance from mid-cell on length axis (", micron(), ")")),
      !have_lenght & have_width  ~ ggplot2::xlab(paste("distance from mid-cell on width axis (", micron(), ")"))
    )
  }

  return(outplot)
}