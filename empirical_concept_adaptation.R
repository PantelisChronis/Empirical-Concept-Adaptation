#Algorithm ECA can be invoked using the function eca_predict()
#Algorithm ADWIN can be invoked using the function adwin_predict()
#The experiments on the synthetic dataset can be performed using the function run_experiments_synth()
#The experiments on the electricity dataset can be performed using the function run_experiments_elec()
#The features for the electricity dataset, required for invoking run_experiments_elec(), are generated using function make_elec_features()

require(cpm)
require(glmnet)
require(e1071)

filter_constant_columns = function(X1,X2){
  X1f = c()
  X2f = c()
  for(i in 1:ncol(X1))
    if(length(unique(X1[,i]))>1){
      X1f = cbind(X1f,X1[,i])
      X2f = cbind(X2f,X2[,i])
    }
  return(list(X1f,X2f))
}

train_val = function(X,y,train_ind,test_ind,type){
  
  mses = c()
  
  X_train = X[train_ind,]
  y_train = y[train_ind]
  
  X_test = X[test_ind,]
  y_test = y[test_ind]
  
  #lam  = lambda_cv(x[train_start:train_end,],y[train_start:train_end])
  if(length(unique(y_train))>1){
  
    Xf = filter_constant_columns(X_train,X_test)
    X_train = Xf[[1]]
    X_test = Xf[[2]]
    
    if(type == 'enet'){ 
      mdl  = glmnet(X_train,y_train,lambda=c(10),alpha=0.01)
      se = (y_test-predict(mdl,X_test))**2
    }
    if(type == 'svm'){
      mdl = svm(X_train,y_train,tolerance = 0.1,kernel='polynomial')
      se = (y_test-predict(mdl,X_test))**2
    }
   
  }
  else{

    y_pred = y_train[1]
    se = (y_test-y_pred)**2
    
  }
  
  return(se)

}

train_val_resamp = function(X,y,train_start,val_end,val_window_size,val_size,type){
 
  ers = c()  
  n = length(y)
  val_win_start = val_end-val_window_size+1
  train_end = val_win_start-1
  all_indices = train_start:val_end
  val_intr_indices = t(matrix(val_win_start:val_end,nrow=val_size))-train_start+1
  times = val_window_size/val_size
  all_test_inds= c()
  
  for(i in 1:times){
    train_ind = all_indices[-val_intr_indices[i,]]
    test_ind = all_indices[val_intr_indices[i,]]
  # all_test_inds = c(all_test_inds,test_ind)
    ers = c(ers, train_val(X,y,train_ind,test_ind,type))
  }
  
  return(ers)
  
}

seq_train_val = function(X,y,curr_time,min_train_size,val_window,val_size,step,type){

  ers = c()
  starting_points = seq(curr_time-min_train_size+1,1,-step)
  
  for(start_point in starting_points)
    ers = c(ers,train_val_resamp(X,y,start_point,curr_time,val_window,val_size,type) )
  
  ers = t(matrix(ers,nrow=val_window))
  plot(starting_points,rowMeans(ers))
  
  return(list(ers,starting_points))
  
}

select_start_from_seq = function(ers,starting_points,alpha){
 
  minInd = which.min(rowMeans(ers))
  pval = wilcox.test(ers[minInd,],ers[nrow(ers),])$p.value
  if (is.na(pval))
    pval=1
  if(pval<alpha)
    return(starting_points[minInd])
  else
    return(1)
  
}

eca_predict = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha,type){
  
  seq = seq_train_val(X,y,curr_time,min_train_size,val_window,val_size,step,type)
  train_start = select_start_from_seq(seq[[1]],seq[[2]],alpha)
  plot(seq[[2]],sqrt(rowMeans(seq[[1]])))
  abline(v=train_start,col=2)
#  print(train_start)
  test_ind = (curr_time+1):(curr_time+test_size)
  train_ind = train_start:curr_time
  se = train_val(X,y,train_ind,test_ind,type)
  rmse = sqrt(mean(se))
  return(rmse)

}

adwin_predict = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha,type){
  
  train_start = select_start_adwin(y,curr_time,alpha)
 # print(train_start)
  abline(v=train_start,col=3)
  test_ind = (curr_time+1):(curr_time+test_size)
  train_ind = train_start:curr_time
  se = train_val(X,y,train_ind,test_ind,type)
  rmse = sqrt(mean(se))
  return(rmse)
  
}

select_start_cpd = function(y,curr_time,alpha,test){
  cp = processStream(y[1:curr_time],ARL0=1000,cpmType=test)$changePoints
  if(length(cp)>0)
    return(cp[length(cp)])
  else
    return(1)
}

naive_predict = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,type){
  
  train_start = 1
#  print(train_start)
  abline(v=train_start,col=4)
  test_ind = (curr_time+1):(curr_time+test_size)
  train_ind = train_start:curr_time
  se = train_val(X,y,train_ind,test_ind,type)
  rmse = sqrt(mean(se))
  return(rmse)
  
}

acp_predict = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,type,win){
  
  train_start = curr_time-win
  if(train_start<0)
    train_start=1
 # print(train_start)
  abline(v=train_start,col=5)
  test_ind = (curr_time+1):(curr_time+test_size)
  train_ind = train_start:curr_time
  se = train_val(X,y,train_ind,test_ind,type)
  rmse = sqrt(mean(se))
  return(rmse)
  
}

compare_methods = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha1,alpha2,type){
  er_eca = eca_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha1,type)
  er_adwin = adwin_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha2,type)
  er_naive = naive_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,type)
  return(c(er_eca,er_adwin,0,0,0,0,er_naive))
}

compare_methods_multiple_points = function(X,y,min_train_size,val_window,val_size,step,test_size,alpha1,alpha2,type){
  ers = c()
  norm = c()
  times = seq(min_train_size,length(y)-test_size,test_size)
  #times = c(230,250,300,430,450,500)
  for(i in times){
    ers = c(ers, compare_methods(X,y,i,min_train_size,val_window,val_size,step,test_size,alpha1,alpha2,type))
    norm = c(norm, mean(y[(i+1):(i+test_size)]))
  }
  return(list(matrix(ers,nrow=7),times,norm))
}

run_cpd = function(x,ind=0){
  res = processStream(x,cpmType="Mann-Whitney")
  plot(x,xlab="day",ylab="error")
  abline(v=res$changePoints,col='red')
  abline(v=res$detectionTimes,col='red',lty=2)
  # title(ind)
  return(list(res$changePoints,res$detectionTimes,x))
}

make_dat_vec = function(day){
  day_mat = c()
  for(d in day){
    new_day = 1:7*0
    new_day[d] = 1
    day_mat=rbind(day_mat,new_day)
  }
  return(day_mat)
}

make_day_vec = function(day){
  day_mat = c()
  for(d in day){
    new_day = 1:7*0
    new_day[d] = 1
    day_mat=rbind(day_mat,new_day)
  }
  return(day_mat)
}

make_hour_vec = function(period){
  hour_mat = c()
  for(p in period){
    new_hour = 1:48*0
    new_hour[round(p*47)] = 1
    hour_mat=rbind(hour_mat,new_hour)
  }
  return(hour_mat)
}

make_elec_features = function(elec2,time){
  
  vic_price = matrix(elec2$vicprice,nrow=48)[time,]
  vic_demand = colSums(matrix(elec2$vicdemand[time:length(elec2$vicdemand)],nrow=48))
  nsw_price = matrix(elec2$nswprice,nrow=48)[time,]
  nsw_demand = colSums(matr9.3ix(elec2$nswdemand[time:length(elec2$nswdemand)],nrow=48))
  transf = colSums(matrix(elec2$transfer[time:length(elec2$transfer)],nrow=48))
  date = matrix(elec2$date,nrow=48)[time,]
  day = make_day_vec(matrix(elec2$day,nrow=48)[time,])
  period = make_hour_vec(matrix(elec2$period,nrow=48)[time,])
  
  y = vic_demand[-(1:30)]
  X = cbind(embed(vic_price,30),embed(vic_demand,30),embed(nsw_price,30),embed(nsw_demand,30),embed(transf,30),day[30:nrow(period),],date[30:length(date)])
  
  X = X[-nrow(X),]

  return(list(X,y))
  
}

select_start_adwin = function(y,curr_time,delta){
  y_part = y[1:curr_time]
  y_part = y_part-min(y_part)
  if(max(y_part)!=0)
    y_part = y_part/max(y_part)
  delta_h = delta/length(y_part)
  sig2 = var(y_part)
  for(i in curr_time:2){
    w0=y_part[1:(i-1)]
    w1=y_part[i:curr_time]
    m=1/(1/length(w0)+1/length(w1))
    e_cut = sqrt(2/m*sig2*log(2/delta_h))+2/3/m*log(2/delta_h)
    if(is.na(e_cut))
      e_cut = Inf
    if(abs(mean(w0)-mean(w1))>e_cut)
      return(i)
  }
  return(1)
}

test_BH = function(pvals,q){
  spvals = sort(pvals)
  dec = c()
  m = length(spvals)
  for(i in 1:m){
    if (i < max(which(spvals <= i/m*q)))
      dec = c(dec,1)
    else 
      dec = c(dec,0)
  }
  return(dec)
}

run_all_water = function(users,min_train,val_win,val_size,back_step,test_size,a1,a2,type){
  all_res = c()
  all_norm = c()
  all_times = c()
  for(i in 1:length(users)){
    user = users[[i]]
    X = user[[1]]
    y = user[[2]]
    res = compare_methods_multiple_points(X,y,min_train,val_win,val_size,back_step,test_size,a1,a2,type)
    all_res = c( all_res, res[[1]] )
    all_times = c( all_times, res[[2]] )
    all_norm = c( all_norm, res[[3]] )
   # print(matrix(all_res,nrow=7))
   # print(rowMeans(matrix(all_res,nrow=7)))
    print(paste('end user ',i,sep=''))
    a = matrix(all_res,nrow=7)
    b = all_norm
    print(rowMeans(t(t(a)/abs(b))))
  }
  return(list(matrix(all_res,nrow=7),all_times,all_norm))
}

run_all_synthetic = function(users,min_train,val_win,val_size,back_step,test_size,a1,a2,type){
  all_res = c()
  all_norm = c()
  all_times = c()
  for(i in 1:length(users)){
    user = users[[i]]
    X = user[[1]]
    y = user[[2]]
    res = compare_methods_multiple_points(X,y,min_train,val_win,val_size,back_step,test_size,a1,a2,type)
    all_res = c( all_res, res[[1]] )
    all_times = c( all_times, res[[2]] )
    all_norm = c( all_norm, res[[3]] )
    # print(matrix(all_res,nrow=7))
    # print(rowMeans(matrix(all_res,nrow=7)))
    print(paste('end user ',i,sep=''))
    a = matrix(all_res,nrow=7)
    b = all_norm
    print(rowMeans(t(t(a)/abs(b))))
  }
  return(list(matrix(all_res,nrow=7),all_times,all_norm))
}

run_experiments_synth = function(synth_data){
  
  print('run1')
  all_wat_enet1 = run_all_water(synth_data,30,30,5,5,15,0.01,0.01,'enet')
  print('run2')
  all_wat_enet2 = run_all_water(synth_data,30,30,5,5,15,0.05,0.05,'enet')
  print('run3')
  all_wat_enet3 = run_all_water(synth_data,30,30,5,5,50,0.1,0.1,'enet')
  print('run4')
  all_wat_enet4 = run_all_water(synth_data,30,30,5,5,15,0.2,0.2,'enet')
  print('run5')
  all_wat_enet5 = run_all_water(synth_data,30,30,5,5,15,0.3,0.3,'enet')
  
  print('run6')
  all_wat_svm1  = run_all_water(synth_data,30,30,5,5,15,0.01,0.01,'svm')
  print('run7')
  all_wat_svm2  = run_all_water(synth_data,30,30,5,5,15,0.05,0.05,'svm')
  print('run8')
  all_wat_svm3  = run_all_water(synth_data,30,30,5,5,50,0.1,0.1,'svm')
  print('run9')
  all_wat_svm4  = run_all_water(synth_data,30,30,5,5,15,0.2,0.2,'svm')
  print('run10')
  all_wat_svm5  = run_all_water(synth_data,30,30,5,5,15,0.3,0.3,'svm')
  
  return(list(all_wat_enet1,all_wat_enet2,all_wat_enet3,all_wat_enet4,all_wat_enet5,all_wat_svm1,all_wat_svm2,all_wat_svm3,all_wat_svm4,all_wat_svm5))

}

run_experiments_elec = function(elef){
  
  X = elef[[1]]
  y = elef[[2]]
  
  print('run1')
  elef_enet1 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.01,0.01,'enet')
  print('run2')
  elef_enet2 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.05,0.05,'enet')
  print('run3')
  elef_enet3 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.1,0.1,'enet')
  print('run4')
  elef_enet4 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.2,0.2,'enet')
  print('run5')
  elef_enet5 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.3,0.3,'enet')
  
  print('run6')
  elef_svm1 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.01,0.01,'svm')
  print('run7')
  elef_svm2 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.05,0.05,'svm')
  print('run8')
  elef_svm3 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.1,0.1,'svm')
  print('run9')
  elef_svm4 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.2,0.2,'svm')
  print('run10')
  elef_svm5 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.3,0.3,'svm')
  
  return(list(elef_enet1,elef_enet2,elef_enet3,elef_enet4,elef_enet5,elef_svm1,elef_svm2,elef_svm3,elef_svm4,elef_svm5))
  
}

make_normalized_errors = function(res){
  normalized = list()
  for(i in 1:length(res)){
    if(is.null(res[[i]]))
      normalized[[i]]=NULL
    else{
    norms = colMeans(matrix(res[[i]][[3]],nrow=15))
    part_res = res[[i]][[1]]
    for(j in 1:nrow(part_res))
      part_res[j,]=part_res[j,]/norms
    part_res = part_res[,is.finite(colSums(part_res))]
    normalized[[i]]= part_res
    }
  }
  return(normalized)
}


print_results = function(res){
  
  n = length(res)
  
  for(i in 1:n){
    a = res[[i]]
    print(rowMeans(a[[1]]))
  }
  
  print('-')
  
  for(i in 1:n){
    i=5
    a = res[[i]]
    a = abs(t(t(a[[1]])/a[[3]]))
    a = a[,is.finite(colSums(a))]
    print(rowMeans(a))
  }
  
}
