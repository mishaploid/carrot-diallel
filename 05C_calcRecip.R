# ------------------------------------------------------------------------------
#   Function to calculate reciprocal effects by environment
#   Reference: B. Griffing (1956) Aust. J. Biol. Sci. 9:463-493
#   S. Turner
#   20 February 2017
# ------------------------------------------------------------------------------

# inputs:
#   trait = name of column with phenotype
#   .id = observation # (from mice mids object) (factor)
#   year = name of column with location/year information (factor)
#   cross = column with hybrid (jk = kj) (factor)
#   recip = column with hybrid (jk != kj) (ordered by male and female parent)
#   data = name of dataframe where info is stored

calcRecip <- function(trait, .id, year, cross, recip, data) {
  df <- data.frame(trait = data[, trait], .id = data[, .id], year = data[, year],
                   cross = data[, cross], recip = data[, recip])
  # means for directed crosses
  r_ijbar <- aggregate(trait ~ cross:year, data = df, FUN = mean, na.action = na.pass)
  colnames(r_ijbar) <- c("recip", "year", "r_ji")
  df <- merge(df, r_ijbar, by = c("recip", "year"))
  newdf <- merge(data, df[, c(".id", "r_ji")], by = ".id")
  return(newdf)
}