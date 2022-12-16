# kre: Kernel Ratio Estimation

## Authors
Mao Li, Chang Yan, Ping He

## Reference
Carlos, H. A., Shi, X., Sargent, J., Tanski, S., & Berke, E. M. (2010). Density estimation and adaptive bandwidths: a primer for public health practitioners. International journal of health geographics, 9(1), 1-8.

## Description
Implementation of kernel ratio estimation. Different from kernel density estimation, kernel ratio estimation uses an adaptive bandwidth that adpat to the population background to address the inhomogeneous issue.

## Example code
There are point of interest and population raster file in the text folder (the extent is city of Columbus). Run following codes (change the file path), you will get the kernel ratio estimation

```
r = raster("~/test/coraster.tif")
poi = st_read("~/test/poi/Retail_Facilities.shp")
poi = st_transform(poi,32119)
Output = kernelRatioEstimation(r,poi,12,1,12000,500,'gaussian')

```
