# bunch of examples
# johaGL2022

# useful func ? : 
getsubstr_ending <- function(vectorofstr, substr){
  outv <- c()
  for (i in vectorofstr){
    if (str_ends(i, substr)){
      outv <- c(outv, i)
    }
  }
  return(outv)
}

