kldiv = function(f1,f2)
{
  f1 = as.data.frame(f1)
  f2 = as.data.frame(f2)
  f1$y2 = NA
  f1$j = NA
  j = 1
  n = length(f1$x)
  for(k in 1:length(f1$x))
  {
    while((f2$x[j]<f1$x[k]) & (j<n))
      j = j+1
    f1$y2[k] = f2$y[j]
    f1$j[k] = j
  }
  f1$diff = c(0,diff(f1$x))
  #show(range(f1$diff))
  return(sum(f1$y*f1$diff*log(f1$y/f1$y2)))
}

kldivsym = function(f1,f2)
{
  0.5*(kldiv(f1,f2)+kldiv(f2,f1))
}