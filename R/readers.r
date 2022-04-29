read_oufti <- function(oufti_matfile) {
  oufti_matfile |>
    R.matlab::readMat(oufti_matfile) |>
    extract_oufti_celllist() -> cell_list

    
  return(cell_list)
}

read_oufti_parameters <- function(oufti_matfile) {
  oufti_matfile |>
    R.matlab::readMat(oufti_matfile) -> raw
  raw$p[,,1] |>  unlist()  |> as.list() -> param
  return(param)
}

read_supersegger_cell <- function(supersegger_matfile) {
  supersegger_matfile |> 
    R.matlab::readMat() |>
    extract_supersegger_celllist() -> cell_list

  return(cell_list)
}

read_supersegger_mesh_ <- function(supersegger_matfile) {
  supersegger_matfile |>
    R.matlab::readMat() |>
    extract_supersegger_mesh() -> mesh
  return(mesh)
}

read_supersegger_mesh <- function(supersegger_path) {
  supersegger_path <- fs::path(supersegger_path)
  if (fs::is_dir(supersegger_path)) {
    mesh_files <- fs::dir_ls(supersegger_path, glob = "*.mat")
    purrr::map(mesh_files, \(x) read_supersegger_mesh_(x)) |>
      bind_rows() -> meshes
    return(meshes)
  }
  else if (fs::is_file(supersegger_path) && fs::path_ext(supersegger_path) == "mat") {
    return(read_supersegger_mesh_(supersegger_path))
  }
  stop("Invalid file format")
}

read_morphometrics <- function(morphometrics_matfile){
  morphometrics_list <- R.matlab::readMat(morphometrics_matfile)
  cell_list <- extract_morphometrics(morphometrics_list)
  return(cell_list)
}
