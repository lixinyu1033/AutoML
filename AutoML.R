# Auto Machine Learning Script - Optimized by Xinyu Li on 2025-02-13
# Yusheng Li 课题组创作
# Ensure you have hub_data and group files prepared in advance.
# hub_data: a matrix of preliminary selected features.
# group: a file containing grouping information.
auto_ml_analysis <- function(hub_data, group, output_dir) {
  
  hub_data=hub_data_Insomnia
  group
  
  ###检查是否所有的包都正确安装
  
  install_and_load <- function(packages) {
    # 检查哪些包已安装
    installed_packages <- installed.packages()[, "Package"]
    
    # 需要安装的包
    missing_packages <- packages[!packages %in% installed_packages]
    
    # 如果有缺失的包，安装它们
    if (length(missing_packages) > 0) {
      cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
      install.packages(missing_packages)
    }
    
    # 加载所有包
    sapply(packages, require, character.only = TRUE)
  }
  
  # 需要加载的包
  required_packages <- c(
    "Boruta", "ggplot2", "tibble", "ggsci", "glmnet", "randomForest", 
    "mlbench", "caret", "xgboost", "Matrix", "PRROC", "shapviz", 
    "dplyr", "tidyverse"
  )
  
  # 调用函数检查并加载包
  install_and_load(required_packages)
  
  
  
  
  mycolors <- c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
                         '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF',
                         '#739B57','#EFE685','#446983')
                         
  group <- as.factor(group)
  output_dir="./results/"
  # Check output directory
  if (missing(output_dir) || output_dir == "") {
    stop("Output directory must be specified.")
  }
  set.seed(123)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Validate inputs
  if (missing(hub_data) || missing(group)) {
    stop("hub_data and group must be provided as data frames.")
  }
  
  # Initialize a list to store selected genes
  all_genes <- list()
  
  # Boruta Model
  cat("Running Boruta model...\n")
  Var.Selec<-Boruta(
    group~.,
    hub_data,
    pValue = 0.01, #confidence level. 可信水平
    mcAdj = TRUE, #是否使用Bonferroni调整
    #if set to TRUE, a multiple comparisons adjustment using the Bonferroni method will be applied. 
    maxRuns = 100, #迭代最多次数
    doTrace = 0,#可以选0-3，运行结果显示的详细程度，0不显示轨迹
    holdHistory = TRUE, #如果设置为TRUE，则存储完整的重要性历史记录，并将其作为结果的ImpHistory元素返回。
    getImp = getImpRfZ #用于获取属性重要性的函数。默认值是 getImpRfZ，它从 ranger 包运行随机森林并收集平均降低精度测量的 Z 分数。
  )
  
  pdf(file.path(output_dir, "Model_Boruta_ImpHistory.pdf"), width = 10, height = 8)
  #运行轨迹图
  ##las=0，横坐标标签水平、纵坐标标签垂直
  ##las=1，横坐标标签水平、纵坐标标签水平
  ##las=2，横坐标标签垂直、纵坐标标签水平
  ##las=3，横坐标标签垂直、纵坐标标签垂直
  plotImpHistory(Var.Selec,
                 whichShadow=c(T,T,T),
                 ylab="Z-Scores",
                 las=2)             
  dev.off()
  # Plot Boruta importance history
  pdf(file.path(output_dir, 'Model_Bortua.pdf'), width = 10, height = 8)
  plot(Var.Selec,
       whichShadow=c(F,F,F),#不绘制阴影的特征
       xlab="",
       ylab="Z-Scores",
       #las=1，
       las=1)               
  dev.off()
  #输出确认、未拒绝的等变量     
  getConfirmedFormula(Var.Selec) #获取确认的变量的公式，对应绿色
  getNonRejectedFormula(Var.Selec) #获取未拒绝变量的公式，对应绿色和黄色
  getSelectedAttributes(Var.Selec,withTentative=FALSE) #返回确认的变量
  Bortua_gene <- getSelectedAttributes(Var.Selec,withTentative=FALSE)
  
  write.table(Bortua_gene, file = file.path(output_dir, "Boruta_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$Boruta <- Bortua_gene
  cat("Boruta Selected Genes:\n", paste(Bortua_gene, collapse = ", "), "\n\n")
  
  # LASSO Model
  cat("Running LASSO model...\n")
  x <- as.matrix(hub_data)
  y <- group
  lasso <- glmnet(x, y, family = 'binomial', nlambda = 1000, alpha = 1)
  cvfit <- cv.glmnet(x, y, nfolds = 10, family = "binomial", type.measure = "deviance")
  
  cat(paste0("#最佳lambda值出现","\n\n",cvfit$lambda.min))
  cat(paste0("lambda.1se值","\n\n",cvfit$lambda.1se))
  coefficients<-coef(cvfit,s=cvfit$lambda.min)#最佳时筛选基因
  Active.Index<-which(coefficients!=0)
  Active.coefficients<- coefficients[Active.Index]
  coefficients
  lasso_gene <- colnames(x)[(Active.Index-1)[-1]]
  lasso_gene 
  #lasso模型的图
  # 提取plot(fit.lasso)中的数据
  #提取lasso中的lambda和回归系数并建立数据框
  plot.data <- data.frame(as.matrix(lasso$lambda), as.matrix(t(lasso$beta)))
  #将款数据格式转换为长数据格式以便作图，更稳定的版本
  plot.data %>%
    {colnames(.)[1] <- "lambda"; .} %>%  # 修改第一列的列名为 "lambda"
    pivot_longer(cols = 2:ncol(.),
                 names_to = "variable",
                 values_to = "coef") -> plot.data
  
  pdf(file.path(output_dir, 'Model_LASSO_coef.pdf'), width = 10, height = 8)
  ggplot(plot.data,aes(log(lambda),coef,color = variable)) + 
    geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
    geom_line(size=1) + 
    xlab("Lambda (log scale)") + 
    #xlab("L1 norm")+
    ylab('Coefficients')+
    theme_bw(base_rect_size = 1.5)+ 
    scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7)))+
    scale_x_continuous(expand = c(0.01,0.01))+
    scale_y_continuous(expand = c(0.01,0.01))+
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=15,color='black'),
          axis.text = element_text(size=12,color='black'),
          legend.title = element_blank(),
          legend.text = element_text(size=12,color='black'),
          legend.position = 'right')+
    annotate('text',x = -5,y=2,label='Optimal Lambda=0.01655567',color='black')+
    guides(col=guide_legend(ncol = 1))# legend cols
  dev.off()
  
  #绘制10x交叉验证的图
  
  cv.df <- data.frame(lambda =cvfit$lambda,#交叉验证中的lambda
                      mse =cvfit$cvm,#mse
                      sd =cvfit$cvsd)#sd
  
  pdf(file.path(output_dir, 'Model_LASSO_cv.pdf'), width = 10, height = 8)
  ggplot(cv.df, aes(log(lambda), mse)) +
    geom_point() + #点图
    geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), width = 0.1) +#添加误差棒
    scale_x_continuous(name = "Log lambda") +
    scale_y_continuous(name = "Mean Squared Error") +
    ggtitle("10-fold Cross-validation using Lasso Regression")+
    geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed", color = "#E41A1C",size=1) +
    geom_vline(xintercept = log(cvfit$lambda.1se), linetype = "dashed", color = "#377EB8",size=1) +
    annotate(geom = "text",label=("lambda.min"),x=-4.5,y=12,color = "#E41A1C",size=3)+
    annotate(geom = "text",label=("lambda.lse"),x=-2.7,y=11.9,color = "#377EB8",size=3)+
    theme_bw()+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  dev.off()
  
  
  
  # Save LASSO genes and print
  lasso_genes <- colnames(x)[which(coef(cvfit, s = cvfit$lambda.min) != 0)[-1]]
  write.table(lasso_genes, file = file.path(output_dir, "LASSO_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$LASSO <- lasso_genes
  cat("LASSO Selected Genes:\n", paste(lasso_genes, collapse = ", "), "\n\n")
  
  
  
  # Random Forest Model
  cat("Running Random Forest model...\n")
  RF <- randomForest(group ~ ., data = hub_data, ntree = 1000, importance = TRUE)
  eror <- importance(RF) %>% as.data.frame()
  imp <- data.frame(ID = rownames(eror), imp = eror$MeanDecreaseAccuracy)
  imp$relative_imp <- (imp$imp - min(imp$imp)) / (max(imp$imp) - min(imp$imp))
  RF_genes <- imp$ID[imp$relative_imp > 0.5]
  rf_data<- imp[order(imp$relative_imp,decreasing = T),c(1,3)]
  
  pdf(file.path(output_dir, 'RF_Model.pdf'), width = 10, height = 8)
  ggplot(rf_data,aes(reorder(ID,relative_imp),relative_imp))+
    geom_bar(stat = "identity",width = 0.2,fill="#1E90FF")+
    geom_point(col="#B22222",size=3)+
    coord_flip()+
    geom_hline(yintercept=0.5,lty=2,col="red",lwd=0.5) +
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.text.x = element_text(colour = "black",size = 9),
          axis.text.y = element_text(colour = "black",size = 9),
          axis.title.x = element_text(colour = "black",size = 12),
          axis.title.y = element_text(colour = "black",size = 12))+
    xlab('Gene')+ylab('Relative importance')
  dev.off()
  
  which.min(RF$err.rate[, 1])
  tmp <- data.frame(trees=1:1000,error=RF$err.rate[,1])
  
  
  pdf(file.path(output_dir, 'Model_RF_trees_OOB.pdf'), width = 10, height = 8)
  ggplot(tmp,aes(trees,error))+
    geom_line(size=0.5,col="gray40")+
    geom_vline(xintercept=35,lty=2,col="red",lwd=0.75) + # 参考线 xintercept = which.min(rf$err.rate[, 1])
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.text.x = element_text(colour = "black",size = 9),
          axis.text.y = element_text(colour = "black",size = 9),
          axis.title.x = element_text(colour = "black",size = 12),
          axis.title.y = element_text(colour = "black",size = 12),
          legend.position = 'top',
          legend.title = element_text(size = 12))+
    xlab('Number of trees')+ylab('OOB error')
  #ggtitle("")
  dev.off()
  
  
  write.table(RF_genes, file = file.path(output_dir, "RF_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$RandomForest <- RF_genes
  cat("Random Forest Selected Genes:\n", paste(RF_genes, collapse = ", "), "\n\n")
  
  
  
  #svm Model
  cat("Svm model...\n")
  group
  control <- rfeControl(functions = caretFuncs,method = "cv", number = 10)  
  
  # 执行SVM-RFE算法
  set.seed(123)
  n <- dim(hub_data)[1]
  results <- rfe(hub_data,            #数据框格式
                 as.factor(group),   #因子格式
                 sizes = c(n:1), 
                 rfeControl = control,
                 method = "svmLinear")
  
  #自动筛选出rfe代码的本质输出就是5个，是由于算法的本质所决定的，下面我将揭露这一本质
  #具体算法见https://github.com/cran/caret/tree/master/R 可以深入理解这一算法的本质内容
  top = 5 #算法默认可修改
  cat("The top ",
      min(top, results$bestSubset),
      " variables (out of ",
      results$bestSubset,
      "):\n   ",
      paste(results$optVariables[1:min(top, results$bestSubset)], collapse = ", "),
      "\n\n",
      sep = "")
  
  pdf(file.path(output_dir, 'SVM_RFE.pdf'), width = 10, height = 8)
  plot(results, type="o",cex.lab = 3,cex.axis = 4)
  dev.off()
  Svm_genes <- results$optVariables[1:min(top, results$bestSubset)]
  
  write.table(Svm_genes, file = file.path(output_dir, "Svm_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$Svm <- Svm_genes
  cat("Svm Selected Genes:\n", paste(Svm_genes, collapse = ", "), "\n\n") 
  
  
  # XGBoost Model
  cat("Running XGBoost model...\n")
  xgboost_hub_data <- hub_data
  as.numeric(table(group)[1])
  as.numeric(table(group)[2])
  
  #将group改为0和1的形式
  xgboost_hub_data$group <- as.numeric(c(rep("0",as.numeric(table(group)[1])),rep("1",as.numeric(table(group)[2]))))
  
  train_matrix <- sparse.model.matrix(group ~ . - 1, data = xgboost_hub_data)
  dtrain <- xgb.DMatrix(data = train_matrix, label = xgboost_hub_data$group)
  res.xgb <- xgboost(data = dtrain, max_depth = 5, eta = 0.5, objective = 'binary:logistic', nround = 25)
  
  # Save XGBoost feature importance and print
  xgb_importance <- xgb.importance(train_matrix@Dimnames[[2]], model = res.xgb)
  xgb_genes <- xgb_importance$Feature
  write.table(xgb_genes, file = file.path(output_dir, "XGBoost_Selected_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$XGBoost <- xgb_genes
  cat("XGBoost Selected Genes:\n", paste(xgb_genes, collapse = ", "), "\n\n")
  
  pdf(file.path(output_dir, 'Model_xgboost.pdf'), width = 10, height = 8)
  ggplot(xgb_importance, aes(x= reorder( Feature,Gain), y=Gain,fill=Feature)) +
    geom_bar(stat="identity") +
    theme_classic() +
    guides(fill=FALSE)+
    #theme(legend.position = )+
    #geom_hline(yintercept = mean(xgb_importance$Gain),lty = 4,col = "darkred",lwd = 0.8) +
    scale_fill_manual(values=mycolors[1:12])+
    coord_flip()+
    theme_bw()+
    ggtitle('XGBoost')+
    theme(plot.title = element_text(size=24,color='black', face = "bold"),
          axis.title.x =element_text(size=18,color='black', face = "bold"),
          axis.text.x =element_text(size=16, color='black', face = "bold"),
          axis.title.y =element_blank(),
          axis.text.y=element_text(size=16,   color='black',face = "bold"),
          legend.title=element_text(size=20, color='black', face = "bold"),
          legend.text=element_text(size=18, color='black', face = "bold"),
          title=element_text(size=20, color='black', face = "bold"),
          strip.text = element_text(size = 14, face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    labs(x="gene",y="Gain",fill="")
  dev.off()
  
  all_genes
  final_gene <- Reduce(intersect, list(
    all_genes$Boruta,
    all_genes$LASSO,
    all_genes$RandomForest,
    all_genes$Svm
  ))
  
  cat("最终筛选的所有模型的交集基因\n","\n\n", final_gene)
  
  
  # Save and print all genes
  all_genes <- as.data.frame(do.call(cbind, lapply(all_genes, function(x) {
    length(x) <- max(sapply(all_genes, length))
    return(x)
  })))
  write.csv(all_genes, file = file.path(output_dir, "all_genes_Genes.csv"))
  # Summary
  cat("Analysis complete. Results are saved in the 'results/' directory.\n")
  #
  
  if (is.null(final_gene) || (length(final_gene) > 0 && is.na(final_gene))) {
    cat("Analysis complete. 'There is no common gene.'\n")
  } else {
    cat("Analysis complete. 'There is no common gene.'\n")
    n1 <- dim(all_genes)[1]
    final_gene <- rep(NA, n1)  # 使用 NA 而非 "NA" 字符串
    all_genes$common_gene <- final_gene
  }
  return(all_genes)
  
}


