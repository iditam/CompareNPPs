library(data.table)
library(doParallel)
library(ROCR)
library(ggplot2)
library(tidyr)
library(cowplot)

num_of_cores <- detectCores()
registerDoParallel(cores=num_of_cores)

ROC_NPP_methods <- function(pos_linkages_file, neg_linkages_file, 
                                 ROC_output_file, PR_output_file, AUCs_output_file, config) {
  # TBD: remove
  # pos_linkages_file = "Linkages/Corum_InComplex_pos_linkages_combinations.csv"
  # neg_linkages_file = "Linkages/Corum_InComplex_neg_linkages_combinations.csv"
  # ROC_output_file = "Results/ROC.pdf"
  # PR_output_file = "Results/PR.pdf"
  
  pos_linkages <- fread(input = pos_linkages_file, verbose = F, showProgress = F)
  neg_linkages <- fread(input = neg_linkages_file, verbose = F, showProgress = F)
  
  if( ("NPPs_from_file" %in% names(config)) && config[["NPPs_from_file"]]) {
    legend_on_bottom <- TRUE
    
    # take details from combinations_details.csv
    combinations_details <- fread(paste0(config$dataFolder, "/", config$combinations_details_file))
    methodsNames <- combinations_details[ , npp_name]
    npp_nice_names <- combinations_details[ , nice_name]
    names(npp_nice_names) <- methodsNames
    npp_colors <- combinations_details[ , color]
    names(npp_colors) <- methodsNames
    
    
  } else {
    legend_on_bottom <- FALSE
    npp_nice_names <- unlist(config$npp_nice_names)
    npp_colors <- unlist(config$npp_colors)
    methodsNames <- config$NPPs
  }

  # deal with multiple correlations:
  # Get all combination of NPP,correlation:
  if (("CORRELATIONs" %in% names(config)) && (length(config[["CORRELATIONs"]]) > 1)) {
    NPP_correlation_comb <- crossing(methodsNames, config[["CORRELATIONs"]])
    colnames(NPP_correlation_comb) <- c("NPP", "correlation")
    NPP_correlation_comb$title = paste0(NPP_correlation_comb$NPP, "-", NPP_correlation_comb$correlation)
    
    # nice_names_comb <- crossing(npp_nice_names[methodsNames], config[["CORRELATIONs"]])
    # colnames(nice_names_comb) <- c("NPP_nice_name", "correlation")
    # 
    # npp_nice_names <- paste0(nice_names_comb$NPP_nice_name, "-", nice_names_comb$correlation)
    # npp_colors <- config[["distinctive_colors"]][1:nrow(nice_names_comb)]
    
    # methodsNames <- NPP_correlation_comb$title
    # npp_nice_names <- sapply(methodsNames, function(mname) {
    #   paste0(npp_nice_names[unlist(NPP_correlation_comb[i,"NPP"])], " & ", 
    #          NPP_correlation_comb[i, "correlation"])
    # })
    # names(npp_nice_names) <- NPP_correlation_comb$title
    # 
    # npp_colors <- sapply(1:nrow(NPP_correlation_comb), function(i) {
    #   npp_colors[unlist(NPP_correlation_comb[i,"NPP"])]
    # })
    # names(npp_colors) <- NPP_correlation_comb$title
    # 
    methodsNames <- NPP_correlation_comb$title
    npp_nice_names <- npp_nice_names[methodsNames]
    names(npp_nice_names) <- methodsNames
    
    npp_colors <- npp_colors[methodsNames]
    names(npp_colors) <- NPP_correlation_comb$title
    
    # cancelled - causes a BUG! as the order of NPP_correlation_comb$title is different.
    # also not needed - as these names already exist from config
    # names(npp_nice_names) <- methodsNames
    # names(npp_colors) <- methodsNames
  } 
  
  allMethodsList <- lapply(methodsNames, function(method) {
    dt <- rbind(pos_linkages[ , .(get(method), Status)], 
                neg_linkages[ , .(get(method), Status)])
    setnames(dt, 1, method)
  })
  names(allMethodsList) <- methodsNames
  

  # calculate ROC
  # auc_vec <- vector("numeric", length(allMethodsList))
  perfAndAUCList <- foreach(i = 1:length(allMethodsList)) %dopar% {
    pred <- prediction(allMethodsList[[i]][ ,1], allMethodsList[[i]][ ,"Status"], label.ordering = c("N", "P"))
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    perfAUC <- performance(pred, measure = "auc")
    auc <- as.numeric(perfAUC@y.values[[1]])
    # perf
    # calculate also Precision-Recall and it's AUC
    PR_perf <- performance(pred, measure = "prec", x.measure = "rec") 
    # x <- PR_perf@x.values[[1]] # Recall values
    # y <- PR_perf@y.values[[1]] # Precision values
    # PR_auc <- sum(diff(x) * (head(y,-1)+tail(y,-1)))/2 # NOT WORKING!
    list(ROC=perf,ROC_auc=auc, PR=PR_perf)
  }
  names(perfAndAUCList) <- names(allMethodsList)
  
  # ROC curve
  rocList <- lapply(perfAndAUCList, function(x) {x[["ROC"]]})
  auc_vec <- sapply(perfAndAUCList, function(x) {x[["ROC_auc"]]})
  # sort the AUCs descending
  auc_vec <- sort(auc_vec, decreasing = T)
  rocDT <- foreach(i = 1:length(rocList), .combine = rbind) %dopar% {
    data.table(FPR=attr(rocList[[i]], "x.values")[[1]],
               TPR=attr(rocList[[i]], "y.values")[[1]],
               Method=rep(names(rocList)[i],
                          length(attr(rocList[[i]], "x.values")[[1]])))
  }
  methods_to_flip <- c("UniProt_no_col_zscore")
  rocDT[Method %in% methods_to_flip, TPR := 2*FPR-TPR]
  auc_vec[methods_to_flip] <- 1-auc_vec[methods_to_flip]
  
  # write AUCs to file
  auc_data <- data.table(method=npp_nice_names[names(auc_vec)],
                         AUC=round(auc_vec, 3))
  fwrite(auc_data, file = AUCs_output_file)
                    
  aucAnno <- paste0(auc_data$method, '=', auc_data$AUC)

  if (legend_on_bottom) {
    pdf(ROC_output_file, height = 10 + length(methodsNames)/2)
  } else {
    pdf(ROC_output_file)
  }
  # png(ROC_output_file, width = 1320, height = 1320, units = "px")
  
  ggROCR_all <- ggplot(data = rocDT, mapping = aes(x = FPR, y = TPR, colour = Method)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, colour="grey", linetype = "longdash") +
    scale_colour_manual(name = 'AUC',
                        labels = aucAnno,
                        breaks = names(auc_vec),
                        values = npp_colors[names(auc_vec)]) +
    labs(x='False positive rate', y='True positive rate',
         title=config$ROC_title) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1")) +
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1")) +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0), 
          legend.margin = margin(5,5,5,5))
  
  if (legend_on_bottom) {
    ggROCR_all <- ggROCR_all +
      theme(legend.position="bottom") +
      guides(colour=guide_legend(ncol = 1,title.position = "top"))
  }
  
  print(ggROCR_all)
  dev.off()
  
  saveRDS(ggROCR_all, file = paste0(ROC_output_file, "_ggplot.rds"))
  # Precision-Recall
  
  PR_List <- lapply(perfAndAUCList, function(x) {x[["PR"]]})
  # auc_vec <- sapply(perfAndAUCList, function(x) {x[["PR_auc"]]})
  # sort the AUCs descending
  # auc_vec <- sort(auc_vec, decreasing = T)
  PR_DT <- foreach(i = 1:length(PR_List), .combine = rbind) %dopar% {
    data.table(X=attr(PR_List[[i]], "x.values")[[1]],
               Y=attr(PR_List[[i]], "y.values")[[1]],
               Method=rep(names(PR_List)[i],
                          length(attr(PR_List[[i]], "x.values")[[1]])))
  }
  
  if (legend_on_bottom) {
    pdf(PR_output_file, height = 13)
  } else {
    pdf(PR_output_file)
  }
  # png(PR_output_file, width = 1320, height = 1320, units = "px")
  gg_PR <- ggplot(data = PR_DT, mapping = aes(x = X, y = Y, colour = Method)) +
    geom_line() +
    scale_colour_manual(name = 'NPP',
                        labels = npp_nice_names[methodsNames],
                        breaks = methodsNames,
                        values = npp_colors[methodsNames]) +
    labs(x='Recall', y='Precision',
         title=config$PR_title
         # , subtitle=paste0("Dataset comprising ", nrow(pos_linkages), 
         #                 " positive gene pairs and ", nrow(neg_linkages), 
         #                 " random negative pairs from ", config$dataset_source)
         ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1")) +
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","0.25","0.5","0.75","1"))
  
  if (legend_on_bottom) {
    gg_PR <- gg_PR +
      theme(legend.position="bottom") +
      guides(colour=guide_legend(ncol = 1,title.position = "top"))
  }

  print(gg_PR)
  dev.off()
  
  saveRDS(gg_PR, file = paste0(PR_output_file, "_ggplot.rds"))
  
}