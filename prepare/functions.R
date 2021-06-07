
# ---- fix-outlier by 99.9 percentile ----
pcap <- function(df, feas){
  x = df[,feas]
  for (i in which(sapply(x, is.numeric))) {
    quantiles <- quantile( x[,i], c(.001, .999 ), na.rm =TRUE)
    x[,i] = ifelse(x[,i] < quantiles[1] , quantiles[1], x[,i])
    x[,i] = ifelse(x[,i] > quantiles[2] , quantiles[2], x[,i])
  }
  df[,feas] = x
  return(df)
}


# ---- subset dataframe by time and columns ----
subset_df <- function(data, safe_id_lst){
  features <- c('age', 'Temp', 'Resp', 'Pulse', 'SBP', 'DBP', 'SpO2',
                'GCS', 'FiO2', 'O2_Flow', 'GLUCOSE', 'TOTAL_BILIRUBIN', 'POTASSIUM',
                'ALBUMIN', 'CALCIUM', 'SODIUM', 'WHITE_BLOOD_CELL_COUNT', 'PHOSPHORUS',
                'PROTIME_INR', 'CREATININE', 'PLATELET_COUNT', 'ALT_GPT', 'CO2',
                'ALKALINE_PHOSPHATASE', 'AST_GOT', 'PCO2', 'CHLORIDE', 'TROPONIN_I',
                'PARTIAL_THROMBOPLASTIN_TIME', 'LACTIC_ACID', 'BLOOD_UREA_NITROGEN',
                'OXYGEN_SATURATION', 'MAGNESIUM', 'FIO2')
  
  # ---- add BSI ----
  if (!'BSI' %in% colnames(data)){ # safe_window datasets don't have BSI labels
    data <- data[which(data$id %in% safe_id_lst), ] # keep safe windows from specified cohort(default pure negative patients)
    data$BSI <- 0
  }
  # ---- add episode ----
  if (! 'episode' %in% colnames(data)){
    data$episode <- "0" # safe window are labeled as episode='0'
  }
  # ---- add relative_time ----
  if (!'relative_time' %in% colnames(data)){
    for (ID in unique(data$id)){
      data[data$id==ID,'relative_time'] = c(1:nrow(data[data$id==ID,]))-97
    }
  }
  # ---- subset dataframe ----
  data <- data %>% select(all_of(features), id, episode, relative_time, BSI, txp) %>% 
    filter(relative_time %in% c(-96:96)) %>%
    as.data.frame()
  
  # ---- add episode id ----
  data$ep_id <- paste(as.character(data$id), as.character(data$episode),sep='_')
  
  # ---- create 4 primary groups ----
  data$group <- "non_neg"
  data$group[which(data$txp==1&data$BSI==1)] <- "txp_pos"
  data$group[which(data$txp==1&data$BSI==0)] <- "txp_neg"
  data$group[which(data$txp==0&data$BSI==1)] <- "non_pos"
  
  # ---- fix factor data type ----
  label_lst <- c('BSI', 'txp', 'group', 'episode')
  data[label_lst] <- lapply(data[label_lst], as.factor)
  
  # ---- add readability for primary episode group----
  data$group.factor <- factor(data$group, levels = c("non_neg", "non_pos", "txp_neg", "txp_pos"))
  levels(data$group.factor) <- c("No-transplant: Negative Cultures and Baseline",
                                 "No-transplant: Positive Cultures",
                                 "Transplant: Negative Cultures and Baseline",
                                 "Transplant: Positive Cultures")
  label(data$group.factor) <- "Primary Groups"
  
  
  # ---- fix-outlier by 99.9 percentile ----
  data <- pcap(data, features)
  
  return(data)
}

# ---- calculate AUC score ----
cstat <- function(p,y=NULL,correct=F) {
  if(class(p)[1]=='lrm' || (length(class(p))>1 & class(p)[2]=='lrm')) {
    y <- p$y
    p <- stats::predict(p)
  }
  if(length(p) != length(y)) stop("Length mismatch")
  j <- !is.na(p) & !is.na(y)
  p <- p[j]; y <- y[j]
  r <- rank(p)
  j <- y==1
  n1 <- as.double(sum(j))
  n2 <- as.double(length(y) - n1)
  R1 <- sum(r[j])
  U <- as.double(R1 - (n1*(n1+1))/2)
  roc <- as.double(U / (n1*n2))
  if(correct && roc<0.5) roc <- 1-roc
  return(roc)
}


# ---- prepare the data frame for modeling ----
# prep_df <- function(imp=NULL, B=-24, A=24, txp_status=NULL){
#   # choose imputation
#   if (is.null(imp)){
#     df = df_original
#   }else if (imp=="mice") {
#       df = df_mice
#   }else if (imp=="median"){
#       df = df_median
#   }
#   # select existing features
#   feas = intersect(colnames(df), predictors)
#   # fix outliers by 99.9 percentile
#   df = pcap(df)
#   # choose time window
#   df = df[df$relative_time > B*4 & df$relative_time < A*4, ]
#   # choose population
#   if(is.null(txp_status)){
#     df_final = df[,c('relative_time','id','ep_id','BSI','txp','episode',feas)]
#   }else if((txp_status==1)|(txp_status==0)){
#     df_final = df[df$txp==txp_status, c('relative_time','id','ep_id','BSI','txp','episode',feas)]
#   }
#   return(df_final)
#   
# }

# ---- balance postive and negative sample size ----
balancing <- function(df, type="under"){
  library(DMwR)
  library(dplyr)
  if (tolower(type)=="smote"){
    df_balanced <- SMOTE(BSI ~ ., df, perc.over = 100, perc.under=200, k = 10)
  }else if (tolower(type)=="over"){
    neg = df[df$BSI==0, ]
    pos = sample_n(df[df$BSI==1,], size=nrow(neg), replace=TRUE)
    df_balanced <- rbind(pos, neg)
  }else{
    pos = df[df$BSI==1, ]
    neg = sample_n(df[df$BSI==0,], size=nrow(pos), replace=FALSE)
    df_balanced <- rbind(pos, neg)
  }
  return(df_balanced)
}

# ---- flatten correlation matrix into dataframe ----
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}

# ---- add topN correlated interaction terms ----
add_x_cols <- function(df, fea1=NULL, fea2=NULL, topN=10){
  if (is.null(fea1) & is.null(fea2)){
    library(ggcorrplot)
    library(Hmisc)
    cor_mtrx <- rcorr(as.matrix(df[,c(predictors)]))
    cor_df <- as.data.frame(flattenCorrMatrix(cor_mtrx$r, cor_mtrx$P))
    fea1 = as.character.factor(cor_df[order(-abs(cor_df$cor))[1:topN],'row'])
    fea2 = as.character.factor(cor_df[order(-abs(cor_df$cor))[1:topN],'column'])
  }
  for (i in 1:min(length(fea1),length(fea2))){
      df$new <- df[,fea1[i]]*df[,fea2[i]]
      names(df)[names(df)=="new"] <- gsub(" ","",paste("x_",fea1[i],"_",fea2[i]))
  }
  print(colnames(df))
  return(df)
}

# ---- bootstrap validation (row-wise) ----
boot_valid_row <- function(df, fml, train_frac=0.8, n=200, blc_type=NULL,pen=0){
  # -- input --
  # df: dataframe name
  # fml: model fml
  # train_frac: fraction training data from the whole dataframe
  # n: times of bootstrap
  # blc_type: type of balancing method
  # pen: penalty term in rms package
  # -- output --
  # res: result matrix including mean/std/se of AIC and AUC
  
  # assert dataframe has column 'id'
  if ('id' %in% colnames(df)){
    
    # balancing data
    if (is.null(blc_type)){
      df_mdl = df
    }else{
      df_mdl = balancing(df, blc_type)
    }
    
    auc = c()
    aic = c()
    for (i in 1:n){
      idx = sample(nrow(df_mdl), nrow(df_mdl)*train_frac)
      df_train = df_mdl[idx,]
      df_val = df_mdl[-idx,]
      mdl = robcov(lrm(fml, data=df_train, x=TRUE,y=TRUE,penalty=pen),cluster = df_train$id)
      y = df_val$BSI
      p = stats::predict(mdl, df_val)
      auc = c(auc, cstat(p,y))
      aic = c(aic, AIC(mdl))
    }
    res = c(pen=pen,round(c(AIC=mean(aic), AIC_sd=sd(aic), AIC_se=qnorm(0.975)*sd(aic)/sqrt(n), AUC=mean(auc), AUC_sd=sd(auc), AUC_se=qnorm(0.975)*sd(auc)/sqrt(n)),4))
    return(res)
  }else{
    return("no patient id")
  }
}

# ---- bootstrap validation (patient-wise) ----
# response must have 2 levels to calculate AUC
boot_valid_patient <- function(df,fml,n=100,n_patient=1,blc_type=NULL,pen=0){
  
  # -- input --
  # df: dataframe name
  # fml: model fml
  # n: times of validation
  # n_patient: number of validating / left-out patients from each group
  # blc_type: type of balancing method
  # pen: penalty term in rms package
  # -- output --
  # res: result matrix including mean/std/se of AIC and AUC
  
  # assert dataframe has column 'id'
  if ('id' %in% colnames(df)){
    
    # create empty aic auc vector
    auc = c()
    aic = c()
    
    # balancing data
    if (is.null(blc_type)){
      df_mdl = df
    }else{
      df_mdl = balancing(df, blc_type)
    }
    
    
    # save all pos/neg patient ids in list
    pos_id = unique(df_mdl[df_mdl$BSI==1,]$id)
    neg_id = unique(df_mdl[df_mdl$BSI==0,]$id)
    
    # loop through all patient
    for (i in 1:n){
      
      # get validation patient ID
      val_id = c(sample(pos_id, n_patient),sample(neg_id, n_patient))
      
      # split train val dataset by patient ID
      df_train = df_mdl[!df_mdl$id%in%val_id,]
      df_val = df_mdl[df_mdl$id%in%val_id,]
      
      # train model on training set
      mdl = robcov(lrm(fml, data=df_train, x=TRUE,y=TRUE,penalty=pen),cluster = df_train$ep_id)
      y = df_val$BSI
      p = stats::predict(mdl, df_val)
      auc = c(auc, cstat(p,y))
      aic = c(aic, AIC(mdl))
    }
    res = c(pen=pen,round(c(AIC=mean(aic), AIC_sd=sd(aic), AIC_se=qnorm(0.975)*sd(aic)/sqrt(n), AUC=mean(auc), AUC_sd=sd(auc), AUC_se=qnorm(0.975)*sd(auc)/sqrt(n)),4))
    return(res)
  }else{
    return("no patient id")
  }
}



# ---- N-fold cross-validation (patient-wise) ----
cross_valid <- function(df,fml,response='BSI',N=10,blc_type=NULL,pen=0,by='id',RMS=TRUE){
  # -- input --
  # df: dataframe name
  # fml: model fml
  # N: N-fold
  # blc_type: type of balancing method
  # pen: penalty term in rms package
  # -- output --
  # res: result matrix including mean/std/se of AIC and AUC
  # assert dataframe has column 'id' or 'ep_id'
  if (by %in% colnames(df)){
    
    # create empty aic auc vector
    auc = c()
    aic = c()
    
    # balancing data
    if (is.null(blc_type)){
      df_mdl = df
    }else{
      df_mdl = balancing(df, blc_type)
    }
    
    # save all pos/neg patient ids in list
    pos_id = unique(df_mdl[df_mdl[response]==1,by])
    neg_id = unique(df_mdl[df_mdl[response]==0,by])
    
    # save N folds of ids in both group
    pos_folds = split(pos_id, ceiling(seq_along(pos_id) / (length(pos_id)%/%N)))
    neg_folds = split(neg_id, ceiling(seq_along(neg_id) / (length(neg_id)%/%N)))
    
    for (i in 1:N){
      # get validation patient ID
      val_id = c(unlist(pos_folds[toString(i)]),unlist(neg_folds[toString(i)]))
      
      # split train val dataset by patient ID
      df_train = df_mdl[!df_mdl[,by]%in%val_id,]
      df_val = df_mdl[df_mdl[,by]%in%val_id,]
      
      # train model on training set
      if(RMS){
        mdl = robcov(lrm(fml, data=df_train, x=TRUE,y=TRUE,penalty=pen),cluster = df_train[,by])
      }else{
        mdl = glm(fml,df_train,family=binomial)
      }
      y = df_val[,response]
      p = stats::predict(mdl, df_val)
      auc = c(auc, cstat(p,y))
      aic = c(aic, AIC(mdl))
    }
  
    res = c(pen=pen,round(c(AIC=mean(aic), AIC_sd=sd(aic), AIC_se=qnorm(0.975)*sd(aic)/sqrt(N), AUC=mean(auc), AUC_sd=sd(auc), AUC_se=qnorm(0.975)*sd(auc)/sqrt(N)),4))
    return(res)
  }else{
    if (by =='id'){
      return("no patient id")
    }else if (by=='ep_id'){
      return("no episode id")
    }
  }
}

# violin ggplot DIY function 
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}


diagnose <- function(df,zirz, return_df=FALSE){
  
  # constant from model
  sirs = c(low.Temp=36,up.Temp=38, Pulse=90, low.WHITE_BLOOD_CELL_COUNT=4, up.WHITE_BLOOD_CELL_COUNT=12, Resp=20)
  
  # prediction using temp cutoffs
  df$sirs_t = 0
  df$zirz_t = 0
  df[(df[,'Temp']<sirs['low.Temp'])|(df[,'Temp']>sirs['up.Temp']),'sirs_t']=1
  df[(df[,'Temp']<zirz['low.Temp'])|(df[,'Temp']>zirz['up.Temp']),'zirz_t']=1
  
  # prediction using WBC cutoffs
  df$sirs_w = 0
  df$zirz_w = 0
  df[(df[,'WHITE_BLOOD_CELL_COUNT']<sirs['low.WHITE_BLOOD_CELL_COUNT'])|(df[,'WHITE_BLOOD_CELL_COUNT']>sirs['up.WHITE_BLOOD_CELL_COUNT']),'sirs_w']=1
  df[(df[,'WHITE_BLOOD_CELL_COUNT']<zirz['low.WHITE_BLOOD_CELL_COUNT'])|(df[,'WHITE_BLOOD_CELL_COUNT']>zirz['up.WHITE_BLOOD_CELL_COUNT']),'zirz_w']=1
  
  # prediction using resp cutoffs
  df$sirs_r = 0
  df$zirz_r = 0
  df[(df[,'Resp']>sirs['Resp']),'sirs_r']=1
  df[(df[,'Resp']>zirz['Resp']),'zirz_r']=1
 
  # prediction using Pulse cutoffs
  df$sirs_p = 0
  df$zirz_p = 0
  df[(df[,'Pulse']>sirs['Pulse']),'sirs_p']=1
  df[(df[,'Pulse']>zirz['Pulse']),'zirz_p']=1
  
  df$sirs_BSI = 0
  df$zirz_BSI = 0
  df$sirs_score = rowSums(df[,c('sirs_p','sirs_r','sirs_t','sirs_w')])
  df$zirz_score = rowSums(df[,c('zirz_p','zirz_r','zirz_t','zirz_w')])
  df[df$sirs_score>=2,'sirs_BSI']=1
  df[df$zirz_score>=2,'zirz_BSI']=1
  
  
  TP = sum(df$BSI==1&df$sirs_BSI==1)
  FP = sum(df$BSI==0&df$sirs_BSI==1)
  TN = sum(df$BSI==0&df$sirs_BSI==0)
  FN = sum(df$BSI==1&df$sirs_BSI==0)
  sensitivity = TP/(TP+FN)
  specificity = TN/(TN+FP)
  ppv = TP/(TP+FP)
  npv = TN/(TN+FN)
  fpr = FP/(FP+TN)
  acc = (TP+TN)/(TP+TN+FP+FN)
  f1 = 2*sensitivity*ppv/(sensitivity+ppv)
  auc_cut = cstat(df$sirs_BSI,df$BSI)
  auc_score = cstat(df$sirs_score,df$BSI)
  res_sirs = c(TP, FP, TN, FN, sensitivity, specificity, ppv, npv, fpr, acc, f1, auc_cut, auc_score)
  
  TP = sum(df$BSI==1&df$zirz_BSI==1)
  FP = sum(df$BSI==0&df$zirz_BSI==1)
  TN = sum(df$BSI==0&df$zirz_BSI==0)
  FN = sum(df$BSI==1&df$zirz_BSI==0)
  sensitivity = TP/(TP+FN)
  specificity = TN/(TN+FP)
  ppv = TP/(TP+FP)
  npv = TN/(TN+FN)
  fpr = FP/(FP+TN)
  acc = (TP+TN)/(TP+TN+FP+FN)
  f1 = 2*sensitivity*ppv/(sensitivity+ppv)
  auc_cut = cstat(df$zirz_BSI,df$BSI)
  auc_score = cstat(df$zirz_score,df$BSI)
  res_zirz = c(TP, FP, TN, FN, sensitivity, specificity, ppv, npv, fpr, acc, f1, auc_cut, auc_score)
  
  # combine together
  res = data.frame(rbind(res_sirs, res_zirz))
  rownames(res) <- c('SIRS','ZIRZ')
  colnames(res) <- c('TP', 'FP', 'TN', 'FN', 'sens', 'spec', 'ppv','npv','fpr','acc','f1','auc_cut','auc_score')
  
  if (return_df){
    return(df)
  }else{
    return(res)
  }
}


calculate_zirz <- function(beta_txp,ct_txp,beta_non,ct_non){
  
  sirs = c(low.Temp=36,up.Temp=38, Pulse=90, low.WHITE_BLOOD_CELL_COUNT=4, up.WHITE_BLOOD_CELL_COUNT=12, Resp=20)
  fea = 'Temp'
  y1 = beta_non[fea]*(sirs['low.Temp']-ct_non[fea])^2
  y2 = beta_non[fea]*(sirs['up.Temp']-ct_non[fea])^2
  zirz = c(low=round(-sqrt( y1/beta_txp[fea] )+ct_txp[fea],1),up=round(sqrt( y2/beta_txp[fea] )+ct_txp[fea],1))
  if(is.na(zirz['low.Temp']) | is.na(zirz['up.Temp']) ){
    zirz['low.Temp'] = sirs['low.Temp'] 
    zirz['up.Temp'] = sirs['up.Temp'] 
  }
  
  fea = 'Pulse'
  y = beta_non[fea]*(sirs[fea]-ct_non[fea])
  zirz = c(zirz,round(y/beta_txp[fea]+ct_txp[fea]),0)
  if(is.na(zirz[fea])){
    zirz[fea] = sirs[fea]
  }
  
  fea = 'WHITE_BLOOD_CELL_COUNT'
  y1 = beta_non[fea]*(sirs['low.WHITE_BLOOD_CELL_COUNT']-ct_non[fea])^2
  y2 = beta_non[fea]*(sirs['up.WHITE_BLOOD_CELL_COUNT']-ct_non[fea])^2
  zirz = c(zirz,c(low=round(-sqrt( y1/beta_txp[fea] )+ct_txp[fea],0), up=round(sqrt( y2/beta_txp[fea] )+ct_txp[fea]),0))
  if(is.na(zirz['low.WHITE_BLOOD_CELL_COUNT']) | is.na(zirz['up.WHITE_BLOOD_CELL_COUNT']) ){
    zirz['low.WHITE_BLOOD_CELL_COUNT'] = sirs['low.WHITE_BLOOD_CELL_COUNT'] 
    zirz['up.WHITE_BLOOD_CELL_COUNT'] = sirs['up.WHITE_BLOOD_CELL_COUNT'] 
  }
  
  fea = 'Resp'
  y = beta_non[fea]*(sirs[fea]-ct_non[fea])
  zirz = c(zirz,round(y/beta_txp[fea]+ct_txp[fea],0))
  if(is.na(zirz[fea])){
    zirz[fea] = sirs[fea]
  }
  
  
  return(zirz)
}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

