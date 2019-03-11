# Calculate runoff based on a (slope adjusted) CN raster from run-off
calc_runoff <- function(Sret, P){
  return( (P > (0.2 * Sret)) * ((P - 0.2 * Sret)^2 / (P + 0.8 * Sret)))
}

# Calculate runoff based on a (slope adjusted) CN raster from run-off
calc_runoff2 <- function(Sret, P){
  return((P - 0.2 * Sret)^2 / (P + 0.8 * Sret))
}

calc_leaching <- function(P){
  return((0.0789 * P) - 0.3594)
}