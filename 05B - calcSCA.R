# ------------------------------------------------------------------------------
#   Function to calculate specific combining ability (SCA) by environment
#   Reference: B. Griffing (1956) Aust. J. Biol. Sci. 9:463-493
#   S. Turner
#   20 February 2017
# ------------------------------------------------------------------------------

calcSCA <- function(trait, .id, year, female, male, cross, data) {
  df <- data.frame(trait = data[, trait], .id = data[, .id], year = data[, year],
                   female = data[, female],
                   male = data[, male], cross = data[, cross])
  # grand mean
  X..bar <- aggregate(trait ~ year, data = df, FUN = mean)
  colnames(X..bar) <- c("year", "X..bar")
  # mean for parent i (as female)
  X_i. <- aggregate(trait ~ female:year, data = df, FUN = mean)
  colnames(X_i.) <- c("female", "year", "X_i.")
  # mean for parent j (as male)
  X_.j <- aggregate(trait ~ male:year, data = df, FUN = mean)
  colnames(X_.j) <- c("male", "year", "X_.j")
  # mean of parent i (as male)
  X_.i <- aggregate(trait ~ female:year, data = df, FUN = mean)
  colnames(X_.i) <- c("male", "year", "X_.i")
  # mean of parent j (as female)
  X_j. <- aggregate(trait ~ male:year, data = df, FUN = mean)
  colnames(X_j.) <- c("female", "year", "X_j.")
  df <- merge(df, X..bar, by = "year")
  df <- merge(df, X_i., by = c("female", "year"))
  df <- merge(df, X_.j, by = c("male", "year"))
  df <- merge(df, X_.i, by = c("male", "year"))
  df <- merge(df, X_j., by = c("female", "year"))
  df$s_ij <- df$trait - 1/2*(df$X_i. + df$X_.j + df$X_.i + df$X_j.) + df$X..bar
  newdf <- merge(data, df[, c(".id", "s_ij")], by = ".id")
  return(newdf)
}
