extract_oufti_celllist <- function(oufti_list) {
  tibble(oufti_list$cellList[[2]]) |>
    unnest(everything()) |>
    rename(cell = 1) |>
    rowid_to_column("frame") |>
    rowwise() |>
    mutate(cell = list(t(cell))) -> cellid

  tibble(oufti_list$cellList[[1]]) |>
    unnest(everything()) |>
    rename(celldata = 1) |>
    rowid_to_column("frame") |>
    inner_join(cellid, by = c("frame" = "frame")) |>
    unnest_longer(c(celldata, cell)) |>
    unnest(celldata) |>
    unnest_wider(celldata) |>
    mutate(
      across(where(\(x) length(dim(x)) == 2 && dim(x)[2] == 1), \(x) c(x)),
      box = st_box(bot = box[,2], top = box[,4]+ box[,2], left = box[,1], right = box[,1]+box[,3])) -> cell_list
    
  cell_list |>
    select(frame, cell, mesh) |> 
    unnest_cell_meshes(mesh) |> 
    get_cell_max_length_and_width(x0, y0, x1, y1, c(cell, frame)) |> 
    pivot_xy_longer(c(x0, x1, y0, y1), num, n) |> 
    summarise(
      mesh = list(st_polygon_autoclose(x, y)),
      max.length = mean(max.length),
      max.width= mean(max.width)) -> meshes
    inner_join(cell_list, meshes, by = c("frame" = "frame", "cell" = "cell")) -> cell_list

  return(celllist)
}

extract_supersegger_celllist <- function(supersegger_list) {
  lookup <- c('cell_id', 'cell_birth_time', 'cell_death_time', 'cell_age', 'fluor1_sum', 'fluor1_mean', 'fluor1_sum_death', 'fluor1_mean_death', 'fluor2_sum', 'fluor2_mean', 'fluor2_sum_death', 'fluor2_mean_death', 'mother_id', 'daughter1_id', 'daughter2_id')
  names(lookup) <- c("cell", "birth", "death", "edgelength", "fluorsum", "fluormean", "fluorsum_D", "fluormean_D", "fluorsum2", "fluormean2", "fluorsum_D2", "fluormean_D2", "parent", "child1", "child2")
  
  as_tibble(supersegger_list$data) -> cell_list
  colnames(cell_list) <- janitor::make_clean_names(unlist(supersegger_list$def))
  cell_list |>
    select(any_of(lookup)) |>
    rename(any_of(lookup)) |>
    mutate(
      parent = if_else(is.na(parent), 0, parent)
    ) -> cell_list
  
}

read_oufti <- function(oufti_matfile){
  oufti_matfile |>
    R.matlab::readMat(oufti_matfile) |>
    extract_oufti_celllist() -> cell_list

    
  return(cell_list)
}

read_supersegger <- function(supersegger_matfile) {
  supersegger_matfile |> 
    R.matlab::readMat() |>
    extract_supersegger_celllist() -> cell_list
  
  return(cell_list)
}


extract_mesh <- function(celllist) {
  celllist |>
    unnest_cell_meshes(mesh, cell, frame) |>
    get_cell_max_length_and_width(x0, x1, y0, y1, c(cell, frame)) |>
    pivot_xy_longer(c(x0, y0, x1, y1), num, n) |>
  return(meshes)
}


unnest_cell_meshes <- function(celllist, mesh) {
  celllist |>
    rowwise() |>
    filter(length({{mesh}}) >= 4) |>
    mutate(mesh = list(tidy_oufti_mesh({{mesh}}))) |>
    unnest({{mesh}}) -> meshes
  return(meshes)
}

tidy_oufti_mesh <- function(mesh) {
  as_tibble(mesh) |>
    rename(x0 = 1, y0 = 2, x1 = 3, y1 = 4) |>
    mutate(
        num = 0:(nrow(mesh) - 1)
      ) -> mesh
  return(mesh)
}

get_cell_max_length_and_width <- function(meshlist, x0, y0, x1, y1, .group) {
  meshlist |>
    group_by(across({{.group}})) |>
    mutate(
      xdist = {{x0}} - {{x1}},
      ydist = {{y0}} - {{y1}},
      widths = polar_distance(xdist, ydist),
      max.width = max(widths),
      xdistL0 = {{x0}} - lag({{x0}}, default = first({{x0}})),
      ydistL0 = {{y0}} - lag({{y0}}, default = first({{y0}})),
      distL0 = polar_distance(xdistL0, ydistL0),
      angL0 = polar_angle(ydistL0, xdistL0),
      rn = row_number(),
      n = n(),
      angd = if_else(
        rn == n,
        polar_angle(lag(ydist), lag(xdist)),
        polar_angle(ydist, xdist))
      ) |> 
    rowwise() |>
    mutate(
      anglength = if_else(
        rn == n,
        angle_from_vertical(abs(angL0), abs(angd)),
        angle_from_horizontal(abs(angL0), abs(angd))),
      steplength = if_else(
        rn == 1,
        0,
        sin(anglength) * distL0)) |>
    group_by(across({{.group}}))  |>
    mutate(
      length = cumsum(steplength),
      max.length = sum(steplength)) |>
      select(-c(xdist, ydist, xdistL0, ydistL0, distL0, angL0, angd, anglength, steplength, length)) -> meshlist
  return(meshlist)
}


add_mesh_pixel_to_unit_conversion <- function(meshlist, pixel_to_unit_ratio, unit) {
  meshlist |>
    mutate(
      max_length_unit = max.length * pixel_to_unit_ratio,
      max_width_unit = max.width * pixel_to_unit_ratio,
      x_rotated_unit = x_rotated * pixel_to_unit_ratio,
      y_rotated_unit = y_rotated * pixel_to_unit_ratio,
      ratio = pixel_to_unit_ratio,
      unit = unit) -> meshlist
    return(meshlist)
}

polygonize_and_get_centroids <- function(meshlist, x, y, .group) {
  meshlist |>
    group_by(across({{.group}}))  |>
    summarise(
      cell_polygon = list(st_polygon_autoclose(x = {{x}}, y = {{y}})),
      bounding_box_angle_centroid = list(get_minimum_bounding_box_centroid_and_angle({{x}}, {{y}}))) |>
    rowwise() |>
    mutate(
      bounding_box_centroid = list(bounding_box_angle_centroid[["centroid"]]),
      bounding_box_angle = (180 - bounding_box_angle_centroid[["angle"]]) * pi / 180,
      angle = if_else(is_polygon_rotation_y_positive(cell_polygon, bounding_box_angle, bounding_box_centroid), bounding_box_angle + pi, bounding_box_angle),
      cell_centroid = list(sf::st_centroid(cell_polygon)),
      area = sf::st_area(cell_polygon)) |>
    select(-c(bounding_box_angle_centroid)) -> meshlist
  return(meshlist)
}

polygonize_and_get_centroids <- function(meshlist, x, y, .group) {
  meshlist |>
    group_by(across({{.group}}))  |>
    summarise(
      cell_polygon = list(st_polygon_autoclose(x = {{x}}, y = {{y}})),
      bounding_box_angle_centroid = list(get_minimum_bounding_box_centroid_and_angle({{x}}, {{y}}))) |>
    rowwise() |>
    mutate(
      bounding_box_centroid = list(bounding_box_angle_centroid[["centroid"]]),
      bounding_box_angle = (180 - bounding_box_angle_centroid[["angle"]]) * pi / 180,
      angle = if_else(is_polygon_rotation_y_positive(cell_polygon, bounding_box_angle, bounding_box_centroid), bounding_box_angle + pi, bounding_box_angle),
      cell_centroid = list(sf::st_centroid(cell_polygon)),
      area = sf::st_area(cell_polygon)) |>
    select(-c(bounding_box_angle_centroid)) -> meshlist
  return(meshlist)
}


rotate_polygons <- function(meshlist, polygon, angle, around) {
  meshlist |>
    mutate(
      rotated_polygon = list(st_rotate_around(polygon = {{polygon}}, theta = {{angle}}, around = {{around}}))) -> meshlist
  return(meshlist)
}
      
    

center_and_rotate_meshes <- function(meshlist, x, y, .group) {
  meshlist |>
    group_by(across({{.group}}))  |>
    mutate(cell_polygon = list(st_polygon_autoclose(x = {{x}}, y = {{y}})),
           bounding_box_angle_centroid = list(get_minimum_bounding_box_centroid_and_angle({{x}}, {{y}}))) |>
    rowwise() |>
    mutate(
           cell_centroid = list(sf::st_centroid(cell_polygon)),
           area = sf::st_area(cell_polygon),
           bounding_box_centroid = list(bounding_box_angle_centroid[["centroid"]]),
           angle = (180 - bounding_box_angle_centroid[["angle"]]) * pi / 180,
           x_mid = bounding_box_centroid[1],
           y_mid = bounding_box_centroid[2],<
           x_centered = {{x}} - x_mid,
           y_centered = {{y}} - y_mid) |>
    group_by(across({{.group}}))  |>
    mutate(
           is_positive_y_rotation = if_else(rotate(x_centered, y_centered, angle)[["y"]][[1]] > 0, TRUE, FALSE),
           angle = if_else(
             is_positive_y_rotation,
             angle + pi,
             angle),
           x_rotated = rotate(x_centered, y_centered, angle)[["x"]],
           y_rotated = rotate(x_centered, y_centered, angle)[["y"]]) |>
           select(-c(is_positive_y_rotation, cell_polygon, bounding_box_angle_centroid)) -> meshlist
  return(meshlist)
}

pivot_xy_longer <- function(celllist, .xy_cols, num, n) {
    celllist |> 
    pivot_longer({{.xy_cols}}, names_to = c(".value", "XY"), names_pattern = "(.)(.)") |>
    mutate({{num}} := if_else(XY == 1, {{n}} - {{num}}, {{num}})) |>
    arrange(XY, {{num}}) |>
    select(-c(XY))-> celllist
    return(celllist)
}

compute_signal_mean_and_sd <- function(celllist) {
  celllist |>
    rowwise() |>
    mutate(
      across(starts_with("signal"), \(x) if_else(is.numeric(x), list(x), list(NA)), .names = "{.col}")) |>
    mutate(
      across(starts_with("signal"), \(x) mean(x, na.rm = TRUE), .names = "mean.{.col}"),
      across(starts_with("signal"), \(x) sd(x, na.rm = TRUE), .names = "sd.{.col}")) -> celllist

  return(celllist)
  }


trim_orphan <- function(datasegger) {
  datasegger |>
    rowwise() |>
    filter(any(parent != 0, !is.na(child1), !is.na(child2))) -> datasegger

  return(datasegger)
}
