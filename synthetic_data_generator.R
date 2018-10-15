many_synthetic1 = function(num=100){
  synth = list()
  for(i in 1:num)
    synth[[i]] = single_synthetic1()
  return(synth)
}

many_synthetic2 = function(num=100){
  synth = list()
  for(i in 1:num)
    synth[[i]] = single_synthetic2()
  return(synth)
}

single_synthetic1 = function(){
  
  x = matrix(runif(100*600)*10,ncol=100)
  y = 1:600*0
  
  w  = runif(100)*10
  w0 = w + runif(100)*3
  w1 = w + runif(100)*0.3
  
  y[1:200]    = x[1:200,]%*%w0+rnorm(200)*100
  y[201:400]  = x[201:400,]%*%w+rnorm(200)*100
  y[401:600]  = x[401:600,]%*%w1+rnorm(200)*100
  
  return(list(x,y))
  
}


single_synthetic2 = function(){
  
  n = 400
  m = 100
  
  x  = matrix(runif(m*n)*10,ncol=m)
  
  w  = runif(100)
  
  y = c()
  ws = c()
  
  for(i in 1:n){
    w = w + rnorm(100)
    ws = rbind(ws, w)
    y = c(y, x[i,]%*%w + rnorm(1)*100 )
  }

  return(list(x,y,ws))
  
}

simul_variance = function(){
  
  
  
}