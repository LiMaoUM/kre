#' Kernel Ratio estimation: adaptive bandwidth kernel density estimation
#' @param pop_grid A raster file contains population info
#' @param point Point of interest data
#' @return A list containing
#' @examples
#' @importFrom raster sf
#' @export
kernelRatioEstimation = function(pop_grid, point, max_bw, min_bw, min_pop)  {
  r <- pop_grid
  p <- point

  # Spatial join: target the raster cell with the point location
  cells <- extract(r,p,cellnumbers=T)

  # Load the cell size of raster
  unit <- res(r)[1]

  # Validate the cell size
  if (unit <= 1){
    warning("The cell seems too small to generate reasonable output, please check the projection")
  }

  # Validate whether the cell is square, if not, show the warning message
  if (unit > res(r)[2]){
    warning("The cell of raster file is not sqare, we will use the shortest length to calculate the bandwidth")
    unit = res(r)[2]
  } else if (unit < res(r)[2]){
    warning("The cell of raster file is not sqare, we will use the shortest length to calculate the bandwidth")
  }

  # Calculate the start and end number of cells for searching
  min_num <- round(min_bw * 1609 / unit) #start from
  max_num <- round(max_bw * 1609 / unit) #end at

  # Create the data.frame for loading the results
  records <- data.frame(col1=c(),col2=c())

  rest_cells <- cells

  # Using for loop searching the bandwidth for each data point

  for (i in min_num:max_num){

    # Number of neighborhood is radius * 2 and plus 1
    num_ber = (i*2)+1

    # Create neighborhood matrix
    neighber <- matrix(1, num_ber, num_ber)
    neighber[(i+1),(i+1)] = 0

    # Filtering adjacent cells for the rest of target data point
    adjacent_cell <- adjacent(r,rest_cells[,"cells"],directions=neighber)
    adjacent_cell <- cbind(adjacent_cell,r[adjacent_cell[,'to']])
    colnames(adjacent_cell)[3] <-  "value"
    row.names(adjacent_cell) <-  adjacent_cell[,"from"]

    # Summing up the population around neighborhood
    cell_sum <-  rowsum(adjacent_cell, row.names(adjacent_cell))

    # Filtering out the cell with enough population support at this level of bandwidth
    met_cells <-  row.names(cell_sum[which(cell_sum[,'value']>min_pop),])

    # Keep the cells without meeting the criteria to next round of filtering
    rest_cells <- rest_cells[which(cell_sum[,'value']<=min_pop),]
    met_cells <- cbind(met_cells,rep(i,length(met_cells)))
    colnames(met_cells) <- c("met_cells", "bw")

    # Record the cells met the criteria
    records <- rbind(records,met_cells)

    # If all cells meet the criteria before reaching out the maximum bandwidth, break the loop
    if (length(rest_cells[,"cells"])<=0){
      break
    }


  }
  return(records)

}
