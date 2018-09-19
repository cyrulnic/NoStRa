
rmse = function(x,y){
  temp <- as.vector(x-y)
  temp <- temp^2
  res <- sqrt(mean(temp))
  return(res)
}