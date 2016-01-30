logMeanExp = function(x){
  biggest_value = max(x)
  log(mean(exp(x - biggest_value))) + biggest_value
}
