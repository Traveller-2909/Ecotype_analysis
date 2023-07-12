#Function to reinterpolate an element to a specified length
interpol <- function(series, length, element){
  ids <- NULL
  values <- NULL
  
  IDOto <- unique(series$ID)
  
  for (i in 1:length(IDOto)) {
    id <- rep(IDOto[i], length.out=length)
    value <- approx(series[which(series$ID==IDOto[i]), element],
                    method = "linear", n=length, rule = 2)[2]
    ids <- c(ids, id)
    values <- unlist(c(values, value))
  }
  
  data.frame(ID=ids, value=values)
  
}