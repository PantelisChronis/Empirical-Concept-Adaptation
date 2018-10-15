#Algorithm ECA can be invoked using the function eca_predict()
#Algorithm ADWIN can be invoked using the function adwin_predict()
#Experiments are run using the run_all_experiments() function

require(cpm)
require(glmnet)
require(e1071)
require(rpart)

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
    
    if(ncol(X_train)<2)
      return(NA)
    
    if(type == 'enet'){ 
      mdl  = glmnet(X_train,y_train,lambda=c(10),alpha=0.01)
      se = (y_test-predict(mdl,X_test))**2
    }
    if(type == 'svm'){
      mdl = svm(X_train,y_train,tolerance = 0.1,kernel='polynomial')
      se = (y_test-predict(mdl,X_test))**2
    }
    if(type == 'tree'){
      df_train = data.frame(y_train,X_train)
      df_test = data.frame(y_test,X_test)
      mdl =  rpart(y_train~.,data=df_train,method="anova")
      #mdl = prune(mdl,mdl$cptable[which.min(mdl$cptable[,"xerror"]),"CP"])
      se = (y_test-as.numeric(predict(mdl,df_test)))**2
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

itm_predict = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha,type){
  
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
  #print(train_start)
  abline(v=train_start,col=4)
  test_ind = (curr_time+1):(curr_time+test_size)
  train_ind = train_start:curr_time
  se = train_val(X,y,train_ind,test_ind,type)
  rmse = sqrt(mean(se))
  return(rmse)
  
}

compare_methods = function(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha1,alpha2,type){
  er_itm = itm_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha1,type)
  er_adwin = adwin_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,alpha2,type)
  er_naive = naive_predict(X,y,curr_time,min_train_size,val_window,val_size,step,test_size,type)
  return(c(er_itm,er_adwin,0,0,0,0,er_naive))
}

compare_methods_multiple_points = function(X,y,min_train_size,val_window,val_size,step,test_size,alpha1,alpha2,type,tms='rand'){
  ers = c()
  norm = c()
  if(tms[1]!="rand")
    times = tms
  else
    times = seq(min_train_size+1,length(y)-test_size, floor((length(y)-test_size-min_train_size-1)/30) )
  #print(times)
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
  for(j in 1:(curr_time-1)){
    w = y_part[j:curr_time]
    delta_h = delta/length(w)
    sig2 = var(w)
    found = F
    for(i in 2:length(w)){
      w0 = w[1:(i-1)]
      w1 = w[i:length(w)]
      m = 1/(1/length(w0)+1/length(w1))
      e_cut = sqrt(2/m*sig2*log(2/delta_h))+2/3/m*log(2/delta_h)
      if(is.na(e_cut))
        e_cut = Inf
      if(abs(mean(w0)-mean(w1))>e_cut)
        found = T
    }
    if(found == F)
      break
  }
  return(j)
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
    #print(matrix(all_res,nrow=7))
    #print(rowMeans(matrix(all_res,nrow=7)))
    print(paste('end user ',i,sep=''))
    a = matrix(all_res,nrow=7)
    b = all_norm
    print(rowMeans(t(t(a)/abs(b))))
  }
  return(list(matrix(all_res,nrow=7),all_times,all_norm))
}

run_experiments_wat = function(good_users){
  
  #  print('run1')
  #  all_wat_enet1 = run_all_water(good_users,30,30,5,5,15,0.01,0.01,'enet')
  #  print('run2')
  #  all_wat_enet2 = run_all_water(good_users,30,30,5,5,15,0.05,0.05,'enet')
  print('run3')
  all_wat_enet3 = run_all_water(good_users,30,30,5,5,50,0.1,0.1,'enet')
  #  print('run4')
  #  all_wat_enet4 = run_all_water(good_users,30,30,5,5,15,0.2,0.2,'enet')
  #  print('run5')
  #  all_wat_enet5 = run_all_water(good_users,30,30,5,5,15,0.3,0.3,'enet')
  
  #  print('run6')
  #  all_wat_svm1  = run_all_water(good_users,30,30,5,5,15,0.01,0.01,'svm')
  #  print('run7')
  #  all_wat_svm2  = run_all_water(good_users,30,30,5,5,15,0.05,0.05,'svm')
  print('run8')
  all_wat_svm3  = run_all_water(good_users,30,30,5,5,50,0.1,0.1,'svm')
  #  print('run9')
  #  all_wat_svm4  = run_all_water(good_users,30,30,5,5,15,0.2,0.2,'svm')
  #  print('run10')
  #  all_wat_svm5  = run_all_water(good_users,30,30,5,5,15,0.3,0.3,'svm')
  
  #  return(list(all_wat_enet1,all_wat_enet2,all_wat_enet3,all_wat_enet4,all_wat_enet5,all_wat_svm1,all_wat_svm2,all_wat_svm3,all_wat_svm4,all_wat_svm5))
  return(list(all_wat_enet3,all_wat_svm3))
  
}

run_experiments_elef = function(elef){
  
  #  print('run4')
  #  all_wat_enet3 = run_all_water(users,30,30,5,5,15,0.3,0.1,'enet')
  
  #  print('run4')
  #  all_wat_svm1 = run_all_water(users,30,30,5,5,15,0.01,0.1,'svm')
  
  X = elef[[1]]
  y = elef[[2]]
  
  print('run1')
  elef_enet1 = c()#compare_methods_multiple_points(X,y,30,30,5,5,15,0.01,0.01,'enet')
  print('run2')
  elef_enet2 = c()#compare_methods_multiple_points(X,y,30,30,5,5,15,0.05,0.05,'enet')
  print('run3')
  elef_enet3 = c()#compare_methods_multiple_points(X,y,30,30,5,5,15,0.1,0.1,'enet')
  print('run4')
  elef_enet4 = c()#compare_methods_multiple_points(X,y,30,30,5,5,15,0.2,0.2,'enet')
  print('run5')
  elef_enet5 = c()#compare_methods_multiple_points(X,y,30,30,5,5,15,0.3,0.3,'enet')
  
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

select_good_users = function(users,n){
  
  selected_users = list()
  for(user in users){
    days= user[[2]]
    if(mean(days)>100 & mean(days)<400)
      selected_users[[length(selected_users)+1]]=user
    if(length(selected_users)==n)
      return(selected_users)
  }
  
  return(selected_users)
  
}

extend_water_features = function(users){
  
  ext_users = list()
  for(user in users){
    X=c()
    y=c()
    for(i in 4:length(user[[2]])){
      x = c(user[[1]][(i-3):(i-1),])
      X = rbind(X,x)
      y = c(y,user[[2]][i])
    }
    ext_users[[length(ext_users)+1]] = list()
    ext_users[[length(ext_users)]][[1]] = X
    ext_users[[length(ext_users)]][[2]] = y
  }
  return(ext_users)
  
}

print_results = function(res){
  
  n = length(res)
  out = c()
  
  for(i in 1:n){
    a = res[[i]]
    print(rowMeans(a[[1]]))
    out = rbind(out, rowMeans(a[[1]]))
  }
  
  print('-')
  
  for(i in 1:n){
    a = res[[i]]
    a = abs(t(t(a[[1]])/a[[3]]))
    a = a[,is.finite(colSums(a))]
    print(rowMeans(a))
    out = rbind(out, rowMeans(a))
  }
  return(out)
}

print_results_at = function(res){
  
  times = sort(unique(res[[2]]))
  er = c()
  
  for(t in times){
    sel = res[[2]]==t
    a = res[[1]][,sel]
    b = res[[3]][sel]
    e = abs(t(t(a)/b))
    e = e[,is.finite(colSums(e))]
    if(is.null(dim(e)))
      er = cbind(er,e)
    else
      er = cbind(er,rowMeans(e))
  }
  
  return(list(er,times))
  
}

read_all_dat = function(){
  path1 = '/home/pant/Desktop/model change/data/'
  ys = c()
  labs = c()
  Xs = c()
  for(i in 1:10 ){
    res = read_dat(paste(path1,'batch',i,'.dat',sep=''))
    labs = c(labs,res[[1]])
    ys = c(ys,res[[2]])
    Xs = rbind(Xs,res[[3]])
  }
  return(list(labs,ys,Xs))
}

read_dat = function(path){
  tab = read.table(path,header=F,as.is=T)
  d =matrix( as.numeric(unlist(strsplit(tab[,1],';'))),nrow=2)
  lab = d[1,]
  conc = d[2,]
  X = c()
  for( i in 2:ncol(tab) ){
    d = matrix(as.numeric(unlist(strsplit(tab[,i],':'))),nrow=2)
    X = cbind(X,d[2,])
  }
  return(list(lab,conc,X))
}

run_experiments_sens = function(sens,gas){
  
  X = sens[[3]][sens[[1]]==gas,]
  y = sens[[2]][sens[[1]]==gas]
  
  print('run2')
  elef_enet2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'enet')
  print('run3')
  elef_enet3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'enet')
  print('run4')
  elef_enet4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'enet')
  
  print('run7')
  elef_svm2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'svm')
  print('run8')
  elef_svm3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'svm')
  print('run9')
  elef_svm4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'svm')
  
  print('run12')
  elef_tree2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'tree')
  print('run13')
  elef_tree3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'tree')
  print('run14')
  elef_tree4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'tree')
  
  return(list(elef_enet2,elef_enet3,elef_enet4,elef_svm2,elef_svm3,elef_svm4,elef_tree2,elef_tree3,elef_tree4))
  
}

run_synth = function(synth,min_train,val_win,val_size,back_step,test_size,a1,a2,type,tms){
  all_res = c()
  all_norm = c()
  all_times = c()
  for(i in 1:length(synth)){
    d = synth[[i]]
    X = d[[1]]
    y = d[[2]]
    res = compare_methods_multiple_points(X,y,min_train,val_win,val_size,back_step,test_size,a1,a2,type,tms)
    all_res = c( all_res, res[[1]] )
    all_times = c( all_times, res[[2]] )
    all_norm = c( all_norm, res[[3]] )
    # print(matrix(all_res,nrow=7))
    # print(rowMeans(matrix(all_res,nrow=7)))
    print(paste('end dataset ',i,sep=''))
    a = matrix(all_res,nrow=7)
    b = all_norm
    print(rowMeans(t(t(a)/abs(b))))
  }
  return(list(matrix(all_res,nrow=7),all_times,all_norm))
}

run_experiments_synth = function(synth){
  #synth2:15,synth1:30
  
  print('run2')
  all_wat_enet2 = run_synth(synth,30,30,5,5,30,0.05,0.05,'enet')
  print('run3')
  all_wat_enet3 = run_synth(synth,30,30,5,5,30,0.1,0.1,'enet')
  print('run4')
  all_wat_enet4 = run_synth(synth,30,30,5,5,30,0.2,0.2,'enet')
  
  print('run7')
  all_wat_svm2  = run_synth(synth,30,30,5,5,30,0.05,0.05,'svm')
  print('run8')
  all_wat_svm3  = run_synth(synth,30,30,5,5,30,0.1,0.1,'svm')
  print('run9')
  all_wat_svm4  = run_synth(synth,30,30,5,5,30,0.2,0.2,'svm')
  
  print('run12')
  all_wat_tree2  = run_synth(synth,30,30,5,5,30,0.05,0.05,'tree')
  print('run13')
  all_wat_tree3  = run_synth(synth,30,30,5,5,30,0.1,0.1,'tree')
  print('run14')
  all_wat_tree4  = run_synth(synth,30,30,5,5,30,0.2,0.2,'tree')
  
  return(list(all_wat_enet2,all_wat_enet3,all_wat_enet4,all_wat_svm2,all_wat_svm3,all_wat_svm4,all_wat_tree2,all_wat_tree3,all_wat_tree4))
  
}

run_tree_elef = function(elef){
  
  #  print('run4')
  #  all_wat_enet3 = run_all_water(users,30,30,5,5,15,0.3,0.1,'enet')
  
  #  print('run4')
  #  all_wat_svm1 = run_all_water(users,30,30,5,5,15,0.01,0.1,'svm')
  
  X = elef[[1]]
  y = elef[[2]]
  
  print('run1')
  elef_tree1 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.01,0.01,'tree')
  print('run2')
  elef_tree2 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.05,0.05,'tree')
  print('run3')
  elef_tree3 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.1,0.1,'tree')
  print('run4')
  elef_tree4 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.2,0.2,'tree')
  print('run5')
  elef_tree5 = compare_methods_multiple_points(X,y,30,30,5,5,15,0.3,0.3,'tree')
  
  return(list(elef_tree1,elef_tree2,elef_tree3,elef_tree4,elef_tree5))
  
}

run_experiments_fli = function(fli){
  
  X = fli[,-ncol(fli)]
  y = fli[,ncol(fli)]
  
  print('run2')
  elef_enet2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'enet')
  print('run3')
  elef_enet3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'enet')
  print('run4')
  elef_enet4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'enet')
  
  print('run7')
  elef_svm2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'svm')
  print('run8')
  elef_svm3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'svm')
  print('run9')
  elef_svm4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'svm')
  
  print('run12')
  elef_tree2 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.05,0.05,'tree')
  print('run13')
  elef_tree3 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.1,0.1,'tree')
  print('run14')
  elef_tree4 = compare_methods_multiple_points(X,y,30,30,5,5,30,0.2,0.2,'tree')
  
  return(list(elef_enet2,elef_enet3,elef_enet4,elef_svm2,elef_svm3,elef_svm4,elef_tree2,elef_tree3,elef_tree4))
  
}

mean_diff_pair_mc = function(a,b,r=10000){
  
  n = length(a)
  
  diff0 = a-b
  diffs = c()
  
  for(i in 1:r)
    diffs = c(diffs,mean(diff0 * sample(c(0,1),n,replace = T)))
  
  
  diff0 = mean(a)-mean(b)
  
  pval1 = sum(diffs<diff0)/r
  pval3 = wilcox.test(a-b)$p.val
  pval2 = t.test(a-b)$p.val
  
  hist(diffs)
  abline(v=diff0,lty=2)
  
  return(c(pval1,pval2,pval3))
  
}

calc_acp = function(res){
  
  er1 = res[[1]][3,res[[2]]==230]
  er2 = res[[1]][4,res[[2]]==250]
  er3 = res[[1]][5,res[[2]]==300]
  er4 = res[[1]][3,res[[2]]==430]
  er5 = res[[1]][4,res[[2]]==450]
  er6 = res[[1]][5,res[[2]]==500]
  
  er = c(rbind(er1,er2,er3,er4,er5,er6))
  
  
  norm1 = res[[3]][res[[2]]==230]
  norm2 = res[[3]][res[[2]]==250]
  norm3 = res[[3]][res[[2]]==300]
  norm4 = res[[3]][res[[2]]==430]
  norm5 = res[[3]][res[[2]]==450]
  norm6 = res[[3]][res[[2]]==500]
  
  norm = c(rbind(norm1,norm2,norm3,norm4,norm5,norm6))
  
  return(er/norm)
  
}

run_all_experiments = function(elef_vic,elef_nsw,chem,fli,synth1,synth2,alpha=0.1){
  
  print("elec_vic")
  X = elef_vic[[1]]
  y = elef_vic[[2]]
  
  print("enet")
  elef_vic_enet <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'enet')
  print("svr")
  elef_vic_svr <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'svm')
  print("tree")
  elef_vic_tree <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'tree')
  
  print("elec_nsw")
  X = elef_nsw[[1]]
  y = elef_nsw[[2]]
  
  print("enet")
  elef_nsw_enet <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'enet')
  print("svr")
  elef_nsw_svr <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'svm')
  print("tree")
  elef_nsw_tree <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'tree')
  
  print("chem")
  X = chem[[3]][chem[[1]]==1,]
  y = chem[[2]][chem[[1]]==1]  
  
  print("enet")
  chem_enet <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'enet')
  print("svr")
  chem_svr <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'svm')
  print("tree")
  chem_tree <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'tree')
  
  print("fli")
  X = fli[[1]]
  y = fli[[2]]
  
  print("enet")
  fli_enet <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'enet')
  print("svr")
  fli_svr <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'svm')
  print("tree")
  fli_tree <<- compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'tree')
  
  print("synth1")
  print("enet")
  synth1_enet <<- run_synth(synth1,30,30,5,5,15,alpha,alpha,'enet',c(230,250,300,430,450,500))
  print("svr")
  synth1_svr <<- run_synth(synth1,30,30,5,5,15,alpha,alpha,'svm',c(230,250,300,430,450,500))
  print("tree")
  synth1_tree <<- run_synth(synth1,30,30,5,5,15,alpha,alpha,'tree',c(230,250,300,430,450,500))
  
  print('synth2')
  print("enet")
  synth2_enet <<- run_synth(synth2,30,30,5,5,15,alpha,alpha,'enet',c(200,400))
  print("svr")
  synth2_svr <<- run_synth(synth2,30,30,5,5,15,alpha,alpha,'svm',c(200,400))
  print("tree")
  synth2_tree <<- run_synth(synth2,30,30,5,5,15,alpha,alpha,'tree',c(200,400))
  
  return(list(elef_nsw_enet,elef_nsw_svr,elef_nsw_tree,elef_vic_enet,elef_vic_svr,elef_vic_tree,chem_enet,chem_svr,chem_tree,fli_enet,fli_svr,fli_tree))
  
}

test_a = function(X,y,alphas){
  res = c()
  for(alpha in alphas){
    out = compare_methods_multiple_points(X,y,30,30,5,5,15,alpha,alpha,'enet')
    res = c(res,print_results(list(out))[2,1:2])
  }
  return(matrix(res,nrow=2))
}

