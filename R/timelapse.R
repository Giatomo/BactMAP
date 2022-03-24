# Timelapse averaging

# To get the division percentage/cell and plot the growth per division group

# 18.10.31 Renske van Raaphorst


# division function:


mark_Division <- function(timelapse, division_percentage_cutoff = 0.75) {
  timelapse$division <- 1
  timelapse <- timelapse[order(timelapse$cell, timelapse$frame), ]

  timelapse %>%
    group_by(cell) %>%
    mutate(is_dividing = division_percentage_cutoff * max.length >= lead(max.length)) %>%
    mutate(division = if_else(is_dividing, lag(division), division))

  return(timelapse)
}


# percentage function:

mark_Percentage <- function(timelapse, av = TRUE) {
  division_time <- timelapse$frame %>% 
    group_by(cell, division) %>%
    summarise(
      min_frame = min(frame),
      max_frame = max(frame),
      av_length = mean(max.length),
      length_var = sd(max.length)) %>%
    mutate(division_time = max_frame - min_frame) %>%
    drop_na(length_var)
  
  ## cutoff based on standard deviation of the mean cell length
  cutoff_cell_length <- stats::median(division_time$length_var, na.rm = TRUE) - stats::sd(division_time$length_var, na.rm = TRUE)

  ## cutoff based on the total division time.take out based on stats::sd?
  cutoff_division_time <- stats::median(division_time$division_time)
  # sdtime <- stats::sd(division_time$division_time)

cutoff_division_time %>%
  mutate(
    fulldivision = case_when(
      length_var <= cutoff_cell_length ~ FALSE,
      division_time < 0.5 * cutoff_division_time ~ FALSE,
      division_time > 2 * cutoff_division_time ~ FALSE,
      TRUE ~ TRUE)
  )
  # put the information back into the "test" dataframe
  timelapse <- merge(timelapse, division_time, all = TRUE)
  timelapse$fulldivision[is.na(timelapse$fulldivision)] <- TRUE

  # calculate growth coefficient per cell using stats::lm()
  divL <- unique(timelapse[, c("cell", "division")])
  coeff <- data.frame(
    "cell" = divL$cell,
    "division" = divL$division,
    "coeff" = unlist(lapply(
      c(1:nrow(divL)),
      function(x) stats::lm(timelapse[timelapse$cell == divL$cell[[x]] & timelapse$division == divL$division[[x]], c("max.length", "frame")])[[1]][[2]]
    ))
  )
  medco <- stats::median(coeff$coeff, na.rm = T)

# TODO to complete
# coeff %>% 
#   mutate(
#     growth = case_when(
#       is.na(coeff) ~ TRUE,
#       coeff < 0.5 * stats::median(coeff, na.rm = TRUE) ~ TRUE,
      
#     )
#   )
  
  coeff$growth <- "u"
  coeff$growth[coeff$coeff < 0.5 * medco | is.na(coeff$coeff)] <- "none"
  coeff <- coeff[order(coeff$coeff), ]
  coeff$gn <- c(1:nrow(coeff))
  # not sure what's going on here
  coeff$growth[coeff$growth != "none"] <- cut(coeff$gn[coeff$growth != "none"], breaks = 3)
  coeff$gn <- NULL
  coeff$growth[coeff$growth == "1"] <- "slow"
  coeff$growth[coeff$growth == "2"] <- "med"
  coeff$growth[coeff$growth == "3"] <- "fast"

  timelapse <- merge(timelapse, coeff)

timelapse %>% 
  mutate(percentage = (frame - min_frame) / division_time * 100) %>%
  filter(percentage <= 100)


  # binning: use "cut"
  timelapse$percentage_binned <- cut(timelapse$percentage,
    breaks = 10,
    labels = paste(c(0:9) * 10, c(1:10) * 10, sep = "-")
  )

  timelapse <- timelapse[!is.na(timelapse$cell), ]

  if (av == TRUE) {
    mean_test <- suppressWarnings(stats::aggregate(timelapse[timelapse$fulldivision == TRUE & timelapse$growth != "none", ][colnames(timelapse) != "percentage_binned"],
      by = list("percentage_binned" = timelapse$percentage_binned[timelapse$fulldivision == TRUE & timelapse$growth != "none"]),
      FUN = mean, na.rm = T
    ))
    return(list("timelapse" = timelapse, "mean_by_percentage" = mean_test))
  }

  if (av == FALSE) {
    return(timelapse)
  }
}

# final function for use
#' @export
perc_Division <- function(timelapse, av = TRUE, plotgrowth = TRUE) {
  timelapseS <- unique(timelapse[, c("cell", "frame", "max.length")])
  timelapseS <- mark_Division(timelapseS)
  if (isTRUE(av)) {
    out <- mark_Percentage(timelapseS, av = TRUE)
    out$timelapse <- merge(timelapse, out$timelapse)
  } 
  else {
    out <- list()
    out$timelapse <- mark_Percentage(timelapseS, av = FALSE)
    out$timelapse <- merge(timelapse, out$timelapse)
  }

  if (isTRUE(plotgrowth)) {
    out$plot_growth <- plot_GrowthTime(out$timelapse)
    if (isTRUE(av)) {
      out$plot_avgrowth <- plot_avGrowth(out$timelapse, out$mean_by_percentage)
    }
  }
  return(out)
}


plot_GrowthTime <- function(timelapse, divisionmarked = TRUE, facet_division = TRUE, variable = "max.length") {
  # get division percentage. of course without plotting or we would get caught in a loop!!
  if (divisionmarked == FALSE) {
    timelapse <- perc_Division(timelapse, av = F, plotgrowth = F)$timelapse
  }
  # make column same as column name indicated.
  timelapse$u <- timelapse[, variable]
  timelapse <- timelapse[timelapse$fulldivision == TRUE, ]
  p <- ggplot2::ggplot(
    timelapse[!is.na(timelapse$division), ],
    ggplot2::aes_string(x = "percentage", y = "u", color = "growth")
  ) +
    ggplot2::geom_line(ggplot2::aes_string(group = "cell"), alpha = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::ylab(variable)
  if (facet_division == TRUE) {
    p <- p + ggplot2::facet_wrap(~division)
  }

  return(p)
}


plot_avGrowth <- function(timelapse, mean_by_percentage, variable = "max.length") {
  timelapse$u <- timelapse[, variable]
  mean_by_percentage$u <- mean_by_percentage[, variable]
  mean_by_percentage$percentage_binned_num <- as.numeric(mean_by_percentage$percentage_binned)
  p <- ggplot2::ggplot(timelapse[timelapse$fulldivision == TRUE & timelapse$growth != "none", ], ggplot2::aes_string(x = "percentage_binned", y = "u")) +
    ggplot2::geom_jitter(alpha = 0.8) +
    ggplot2::ylab(variable) +
    ggplot2::geom_line(
      data = mean_by_percentage,
      ggplot2::aes_string(x = "percentage_binned_num", y = "u"), color = "orange", size = 2
    ) +
    ggplot2::theme_classic()
  return(p)
}