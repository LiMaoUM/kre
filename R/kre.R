#' Kernel Ratio estimation: adaptive bandwidth kernel density estimation
#'
#'
#' @param pop_grid A \code{\link[raster]{Raster-class}} contains population info
#' @param point \code{\link[sf]{sf} contains point of interest data
#' @param max_bw maximum searching radius (unit in mile)
#' @param min_bw minimum searching radius (unit in mile)
#' @param min_pop minimum population search radius needs to cover for each point
#' @param cell_size resolution of output raster file (unit in meter)
#'
#'
#' @return \code{\link[raster]{Raster-class}}
#' @export
#'
#' @import raster sf dplyr methods Matrix
#' @importFrom raster extract raster pointDistance coordinates
#' @importFrom sf st_geometry st_union st_convex_hull st_buffer st_as_sf st_coordinates
#' @importFrom terra res
#' @importFrom utils txtProgressBar setTxtProgressBar

kernelRatioEstimation = function(pop_grid,
                                 point,
                                 max_bw,
                                 min_bw,
                                 min_pop,
                                 cell_size,
                                 kernel = c("gaussian","uniform","quartic","triweight","epanechnikov")) {
  r <- pop_grid
  p <- point
  side_offset <- 0 

  # Spatial join: target the raster cell with the point location
  cells <- extract(r,p,cellnumbers=T)

  # Append cell numbers to point data
  p = cbind(p,cells)

  # Load the cell size of raster
  unit <- res(r)[1]

  # Validate the cell size
  if (unit <= 1){

    stop("The cell seems too small to generate reasonable output, please check the projection")

  }

  # Validate whether the cell is square, if not, show the warning message
  if (unit > res(r)[2]){

    warning("The cell of raster file is not sqare, we will use the shortest length to calculate the bandwidth")

    # Update the searching unit
    unit = res(r)[2]

  } else if (unit < res(r)[2]){

    warning("The cell of raster file is not sqare, we will use the shortest length to calculate the bandwidth")

  }

  # Calculate the start and end number of cells for searching
  min_num <- round(min_bw * 1609 / unit) #start from
  max_num <- round(max_bw * 1609 / unit) #end at

  # Set progress bar
  pb  <- txtProgressBar(min_num,
                        max_num,
                        style=3)

  # Create the data.frame for loading the results
  records <- data.frame(col1=c(),
                        col2=c())

  rest_cells <- unique(cells[,'cells'])

  # Using for loop searching the bandwidth for each data point
  message("Start calculating the bandwidth for each point of interest")

  for (i in min_num:max_num){

    setTxtProgressBar(pb, i)

    # Number of neighborhood is radius * 2 and plus 1
    num_ber = (i*2)+1

    # Create neighborhood matrix
    neighber <- matrix(1,
                       num_ber,
                       num_ber)
    neighber[(i+1),(i+1)] = 0

    # Filtering adjacent cells for the rest of target data point
    adjacent_cell <- adjacent(r,
                              rest_cells,
                              directions=neighber)
    adjacent_cell <- cbind(adjacent_cell,
                           r[adjacent_cell[,'to']])
    colnames(adjacent_cell)[3] <-  "value"
    row.names(adjacent_cell) <-  adjacent_cell[,"from"]

    # Summing up the population around neighborhood
    cell_sum <-  rowsum(adjacent_cell,
                        row.names(adjacent_cell))

    # Filtering out the cell with enough population support at this level of bandwidth
    met_cells <-  row.names(cell_sum[which(cell_sum[,'value']>min_pop),,drop=F])

    # Keep the cells without meeting the criteria to next round of filtering
    rest_cells <- as.numeric(row.names(cell_sum[cell_sum[,'value']<=min_pop,,drop=F]))

    # Record the cells met the criteria
    met_cells <- cbind(met_cells,rep(i,length(met_cells)))

    # Name the column
    colnames(met_cells) <- c("met_cells", "bw")

    # Record the cells met the criteria
    records <- rbind(records,met_cells)

    i= i+1

    # If all cells meet the criteria before reaching out the maximum bandwidth, break the loop
    if (length(rest_cells)<=0){
      break
    }


  }

  # Transform the number of cells to meter
  records$bw <- as.numeric(records$bw) * unit

  # Assign the adpative bandwidth to each point of interest
  p <- merge(p,records,
             by.x='cells',by.y='met_cells',all=T)

  # Create the raster file for kernel ratio output
  buffered_poi <- p %>%
    st_geometry() %>%
    st_union() %>%
    st_convex_hull() %>%
    st_buffer(side_offset) %>%
    st_as_sf()

  output_raster <- raster(buffered_poi,
                           resolution = c(cell_size, cell_size))

  # Extract the coordinates from point of interest and output_raster file
  coord_poi <- st_coordinates(p)
  raster_poi <- coordinates(output_raster)

  # Create distance matrix between raster and point of interest (all pairs)
  distmat <- pointDistance(raster_poi,coord_poi,lonlat = F)

  # Force the value larger than distance to 0
  distmat <- ifelse((distmat/p$bw)>=1,0,distmat/p$bw)



  # Convert distance matrix to a sparse matrix
  distmat <- as(distmat, "sparseMatrix")


  # Calculate the kernel ratio
  if (kernel == "gaussian"){

    vec <- GaussianKernelDensity(t(distmat))

  } else if (kernel == "uniform"){

    vec <- UniformKernelDensity(t(distmat))

  } else if (kernel == "quartic") {

    vec <- QuarticKernelDensity(t(distmat))

  } else if (kernel == "triweight"){

    vec <- TriweightKernelDensity(t(distmat))

  } else if (kernel == "epanechnikov"){

    vec <- EpanechnikovKernelDensity(t(distmat))

  }

  # Assign the kernel ratio to each cell of output raster
  output_raster[] <- vec

  return(output_raster)

}
