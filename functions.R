#### Utility Wrappers ####

#' @title Show faults of an expression interactively
#' @description Evaluates the given expression, function and its arguements
#' return the output of the expression. If any errors/warnings are found a message
#' is logged in the console. If used in a shiny App a notificatio is shown, for any
#' errors process from the shiny app is stopped. Warning are shown and the shiny
#' process continues as developed.
#' @param ... Arguements passed to conditions, the expression to execure
#' @param session A shiny object of shiny web app, defaults to NULL
#' @return A return object of type as specified by the expression being executed
#' @references https://github.com/SumedhSankhe/workBench/blob/fd5e808e86e849a933a924735d30348974212fc3/functions.R#L1-L37
show_faults <- function(..., session = NULL){
  warn <- err <- NULL
  res <- withCallingHandlers(tryCatch(..., error = function(e) {
    err <<- conditionMessage(e)
    NULL
  }), warning = function(w) {
    warn <<- append(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  
  if(!is.null(err)){
    if(!is.null(session)){
      shiny::showNotification(as.character(err), duration = 20, type = 'warning',
                              session = session, id = "error") 
      shiny::validate(shiny::need(is.null(err), as.character(err)))
    } else {
      stop(Sys.time(),": ",err)
    }
  } else if (!is.null(warn)){
    warn <- paste(unique(warn), collapse = ", ")
    if(!is.null(session)){
      shiny::showNotification(as.character(warn), duration = 3, type = 'warning',
                              session = session) 
      return(res)
    } else {
      warning(Sys.time(),": ",warn)
      return(res)
    }
  } else {
    return(res)
  }
}


#' @title Display Status messages to the console of the shiny-session
#' @description A utility function which uses the native R and shiny reporting functions
#' to display status of a particular process to the console or the shiny message
#' bar with progress
#' @param detail A character string of message to be displayed inthe console/shiny progress
#' @param value A numeric value 0-1 which depicts the progress of the task at hand
#' @param session A shiny session object
#' @references https://github.com/SumedhSankhe/workBench/blob/fd5e808e86e849a933a924735d30348974212fc3/functions.R#L51
status <- function(detail, value, session = NULL){
  if(!is.null(session))
    shiny::setProgress(value = value, message = "Progress:", detail = detail,
                       session = session)
  message(Sys.time(),": ",detail,"...")
}


#' @title MSstats Theme
#' @description A utility function that standardized all the ggplot themes
#' @param x.axis.size A numeric value for size for the elements on the x-axis
#' @param y.axis.size A numeric value for size for the elements on the y-axis
#' @param legend.size A numeric value fot the size of the legend
theme_MSstats <- function(x.axis.size = 10, y.axis.size = 10, legend.size = 7){
  
  th <- ggplot2::theme(panel.background = element_rect(fill = "white", 
                                                       colour = "black"),
                       panel.grid.major = element_line(colour = "gray95"), 
                       panel.grid.minor = element_blank(), 
                       strip.background = element_rect(fill = "gray95"), 
                       strip.text.x = element_text(colour = c("#00B0F6"), 
                                                   size = 14),
                       axis.text.x = element_text(size = x.axis.size, 
                                                  colour = "black"), 
                       axis.text.y = element_text(size = y.axis.size, 
                                                  colour = "black"),
                       axis.ticks = element_line(colour = "black"), 
                       axis.title.x = element_text(size = x.axis.size + 5,
                                                   vjust = -0.4), 
                       axis.title.y = element_text(size = y.axis.size + 5,
                                                   vjust = 0.3),
                       title = element_text(size = x.axis.size + 4,
                                            vjust = 1.5),
                       legend.key = element_rect(fill = "white", 
                                                 colour = "white"),
                       legend.position = "right", 
                       legend.text = element_text(size = legend.size), 
                       legend.title = element_blank())
  return(th)
}

#### Data Exploration ####



format_summary_table <- function(data = NULL){
  #create crosstable for the conditions vs bioreplicates
  biorep <- unique(data[,.(bioreplicate, condition)])
  biorep <- xtabs(~condition, data = biorep)
  
  #create crosstable for the conditions vs runs if runs data exists
  if('run' %in% names(data)){
    msruns <- unique(data[,.(run, condition)])
    msruns <- xtabs(~condition, data = msruns)
  }else{
    msruns <- rep(0, length(names(biorep)))
    names(msruns) <- names(biorep) #make runs data 0 if not found
  }
  #format it correctly
  summary <- rbind(biorep, msruns)
  rownames(summary) <- c("# of Biological Replicates", "# of MS runs")
  return(summary[,which(colSums(summary, na.rm = T) > 0)])
}

estimate_variance <- function(abundance, annotation){
  
  if(!'data.table' %in% (.packages())){
    library('data.table')
  }
  #check if the data provided is a data.table
  #set all names to lower cases
  if(!'data.table' %in% class(abundance)){
    abundance <- data.table::as.data.table(abundance, keep.rownames = T)
    data.table::setnames(abundance, 'rn', 'protein')
  }
  
  if(any(c('X','V1') %in% names(abundance))){
    data.table::setnames(abundance, c('X', 'V1'), c('protein','protein'),
                         skip_absent = T)
  }
  #check if the data provided is a data.table
  #set all names to lower cases
  if(!'data.table' %in% class(annotation)){
    annotation <- data.table::as.data.table(annotation)
  }
  #standardize names
  names(annotation) <- tolower(names(annotation))
  #convert the data to long form
  abundance <- data.table::melt(abundance, id.vars = 'protein')
  #set the name to identify correctl
  data.table::setnames(abundance, c('variable', 'value'), c('bioreplicate', 'abun'))
  #merge annotation and abundance data
  dt <- merge(annotation, abundance, by = 'bioreplicate')
  dt <- dt[!is.na(abun)]
  groups <- dt[, unique(condition)]
  proteins <- dt[,unique(protein)]
  #split the long data into smaller subsets by proteins
  df_list <- split(dt, by = 'protein')
  #get summary statistics for each protein and condition group
  stat_report <- lapply(df_list, group_summary, group = groups)
  
  return(data.table::rbindlist(stat_report))
}

group_summary <- function(df, group){
  #get the protein names
  prot <- unique(df[,protein])
  #estimate the linear relation between the protein and abundance
  df_lm <- try(lm(abun~condition, data = df), T)
  
  if(!inherits(df_lm, 'try-error')){
    coefs <- df_lm$coefficients #get the coefs of the linear model
    if(length(coefs) == length(group)){
      var <- anova(df_lm)["Residuals", "Mean Sq"] #find the residuals/variance
      var <- rep(sqrt(var), length(group))
      names(var) <- group
      
      coefs[-1] <- sum(coefs, na.rm = T) #calculate the means
      names(coefs) <- gsub('condition','', names(coefs))
      names(coefs)[1] <- setdiff(group, names(coefs))
    }else{
      warning(Sys.time(),": Linear Model coeffienct not equal to n-groups for - ", prot)
      coefs <- rep(NA, length(group))
      names(coefs) <- group
    }
  }else{
    warning(Sys.time(),": Linear Model Failed for protein - ", prot)
    coefs <- rep(NA, length(group))
    names(coefs) <- group
  }
  #format the statistics as a data.table and return
  data.frame('protein' = prot, 'mean_' = t(coefs),
             'sigma_' = t(var), 'sampleMean' = df[,mean(abun, na.rm = T)])
}

meanSDplot <- function (data, smoother_size = 1, xlimUp = 30, ylimUp = 3){
  #plot
  p <- ggplot(data = data, aes(x = mean, y = sigma))+
    stat_density2d(aes(fill = ..density..^0.25), geom = 'tile', contour = F, 
                   n = xlimUp*10)+
    geom_point(alpha = 0.2, shape = 20)+
    scale_fill_continuous(low = "white", high = "#0072B2")+
    geom_line(aes(x = lowess.x, y = lowess.y), color = 'orange',
              size = smoother_size)+
    labs(x = "Mean protein abundance per condition", 
         y = "Standard deviation per condition") +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimUp)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, xlimUp)) +
    theme_MSstats()+
    theme(legend.position = 'none')
  
  return(p)
}


#' @title Formats data to required longer format
#' @description Formats the input data in the required long and wide format
#' depending on input data
#' @param  format A decided which data to use, user provided or default examples
#' @param count A path to the data to be read in, any csv formats are acceptable
#' @param annot An annotation file path as csv
#' @param session A shiny session variable to make the notifications interactive
#' @return A named list of the data an other objects as required
format_data <- function(format, count = NULL, annot = NULL, session = NULL){
  #shiny::validate(shiny::need(format %in% FORMATS, 'Undefined Format'))
  if(format == 'standard'){
    status(detail = 'Importing Protein Abundance file', value = 0.4,
           session = session)
    #read abundance data from the file path provided
    wide <- fread(count$datapath, keepLeadingZeros = T)
    #No column names expected for the protein columns
    #TODO make this name agnostic 
    setnames(wide, 'V1', 'protein')
    uniq_prots <- unique(wide[, protein])
    prots_combinations <- uniq_prots[grep(',|;',uniq_prots)]
    
    if(length(prots_combinations) > 0){
      single_prots <- uniq_prots[-grep(',|;', uniq_prots)]
      v <- do.call('c', lapply(prots_combinations, function(x){
        comb <- unlist(strsplit(x,',|;'))
        val <- any(comb %in% single_prots)
        if(!val){
          return(x)
        }
      }))
      single_prots <- c(single_prots,v)
      message(Sys.time()," Old Proteins", length(uniq_prots),
              ", New Proteins", length(single_prots))
      
      wide <- wide[protein %in% single_prots]
    }
    
    name <- count$name
    status(detail = 'Importing Annotations file', value = 0.5,
           session = session)
    #read annotations from the file path provided
    annot <- fread(annot$datapath, keepLeadingZeros = T)
  }else if(format == 'examples'){
    status(detail = 'Importing Data from MSstatsSampleSize Package', value = 0.5,
           session = session)
    #example data from the package
    wide <- as.data.table(MSstatsSampleSize::OV_SRM_train,
                          keep.rownames = T)
    #examples data from the package
    annot <- as.data.table(MSstatsSampleSize::OV_SRM_train_annotation)
    setnames(wide,'rn','protein')
    name <- "Ovarian Cancer SRM study"
  }else{
    stop('Not Defined')
  }
  
  names(annot) <- tolower(names(annot))
  names(annot)[like(names(annot), 'run|Run|runs|Runs')] <- 'run'
  #get summary table in the right format
  data_summary <- format_summary_table(data = annot)
  
  status(detail = 'Stretching the data', value = 0.6, session = session)
  #convert data to the long form
  data <- melt(wide, id.vars = 'protein', variable.name = 'bioreplicate',
               value.name = 'abundance')
  status(detail = 'Merging Abundance & Spectral Counts with Annotations', 
         value = 0.7, session = session)
  #merge the abundance data with the annotation data with the correct runs/bioreps
  if("run" %in% names(annot)){
    data <- merge(data, annot[, bioreplicate := NULL],
                  by.y = "run", by.x = "bioreplicate")
    #warnings suppressed to avoide it popping up when the column isn't found
    #too lazy to fix this
    suppressWarnings({
      annot[, bioreplicate := NULL]
      data[, bioreplicate.y := NULL]
    })
    setnames(annot, 'run', 'bioreplicate')
  }else{
    data <- merge(data, annot, by = 'bioreplicate') 
  }
  #Get variance and mean estimations for the data
  status(detail = "Estimating Mean & Variances", value = 0.8, session = session)
  var_summary <- estimate_variance(abundance = wide, annotation = annot)
  plot_data <- copy(var_summary)
  
  status(detail = "Formatting Mean Variance Plot Data", value = 0.9, session = session)
  #get rid of unwanted columns
  plot_data[, sampleMean := NULL]
  #convert data to long form
  plot_data <- melt(plot_data, id.vars = 'protein')
  #split the condition and statistics 
  plot_data[, c("stat", "condition") := tstrsplit(variable, "_.", fixed=TRUE)]
  plot_data <- dcast(plot_data, protein + condition ~ stat, value.var = 'value')
  #find the lowess line
  plot_data[, c('lowess.x', 'lowess.y') := lowess(mean, sigma)]
  
  status(detail = "Generating Mean Variance Estimate Plots", value = 0.95,
         session = session)
  mean_sd_plot <- meanSDplot(data = plot_data)
  
  return(list('long_data' = data, 'data_summary' = data_summary,
              'var_summary' = var_summary, 'b_group' = unique(annot[, condition]),
              'n_prot' = nrow(wide[,.N, protein]), 'mean_sd_plot' = mean_sd_plot,
              'n_group' = nrow(annot[,.N,condition]), 'dataset_name' = name))
}


#' @title Do Principal Component Analysis
#' @description A wrapper function to the base `prcomp()` formats the results 
#' in a required output format
#' @param sim_x A dataframe of the feature variables
#' @param sim_y A dataframe of the predictor vairable
#' @return A named list of the outputs 
do_prcomp <- function(sim_x, sim_y){
  result.pca <- prcomp(sim_x, scale. = TRUE)
  summary.pca <- summary(result.pca)
  important.pc <- result.pca$x[, 1:2]
  pc.result <- data.frame(important.pc, group = sim_y)
  exp.var <- summary.pca$importance
  exp.var <- format(exp.var[2, 1:2] * 100, digits = 3)
  return(list("pc.result" = pc.result, "exp.var" = exp.var))
}


#' @title Make pca plots
#' @description  A wrapper function which formats all the data to pass to the 
#' `pca_plot()` function
make_pca_plots <- function(simulations, choice, width = 3, height = 3,
                           session = NULL){
  
  choice <- unlist(strsplit(choice,'Sample'))
  choice[2] <- gsub(' Size ', '', choice[2]) 
  choice[1] <- gsub(' ','', choice[1])
  status(detail = 'Extracting Simulated Data', value = 0.3, session = session)
  val <- sprintf("Sample_%s_Simulation_%s", choice[2], gsub("Simulation",'', choice[1]))
  df <- copy(simulations[[choice[2]]]$simulated[[val]])
  status(detail = 'Performaing PCA on the data', value = 0.7, session = session)
  pr_comp <- do_prcomp(sim_x = df[,-1], sim_y = as.vector(df[,1]))
  p <- pca_plot(data = pr_comp$pc.result, exp_var = pr_comp$exp.var)
  
  return(p)
}

#' @title Plot PCA outputs
#' @description A utility wrapper for plotting pca data as returned by the `do_prcomp()`
#' @param data A data frame containing the Principal components to be plotted
#' @param exp_car A vector containing numeric values of the expected variance
#' @dot_size A aesthetics field to be passed to the ggplot2 functions
pca_plot <- function(data, exp_var, dot_size = 3){
  p <- ggplot(data = data,
              aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = dot_size) + 
    labs(title = "Input dataset",
         x = sprintf("PC1 (%s%% explained var.)",exp_var[1]),
         y = sprintf("PC2 (%s%% explained var.)",exp_var[2])) +
    theme_MSstats()
  return(p)
}



#### Data Simulation ####

#' @title Simulate datasets to be tested out
#' @description A wrapper function for the `simulateDataset` function from the 
#' MSstatsSampleSize package which enables simulating datasets for running experiments
#' @param data 
simulate_grid <- function(data = NULL, stats = NULL, n_group, n_prots, num_simulation,
                          exp_fc, list_diff_proteins = NULL, sel_simulated_proteins ='proportion', 
                          prot_proportion = 1, prot_number = 1000,
                          samples_per_group, sim_valid = F,
                          valid_samples_per_grp = 50, seed = NULL, session = NULL){
  
  status(detail = "Setting Up Data Simulation Runs", value = 0.1, session = session)
  if(!is.null(seed))
    set.seed(seed)
  
  if(exp_fc != 'data'){
    status(detail = "Extracting Fold Change Informations", value = 0.15, session = session)
    diff_prots <- unlist(strsplit(list_diff_proteins, ","))
    fc <- data.table('condition' = exp_fc$orig_group,
                     'exp_fc' = exp_fc$`Fold Change Value`)
  } else{
    diff_prots <- NULL
    fc <- exp_fc
  }
  status(detail = "Extracting Number of Samples Information", value = 0.2, session = session)
  samp <- as.numeric(unlist(strsplit(samples_per_group, ',')))
  shiny::validate(shiny::need(all(!is.na(samp)),
                              sprintf("Samples Per Group need to be numeric values, Found : %s",
                                      samples_per_group)),
                  shiny::need(all(samp >= 1), "All samples Need to be >= 1"))
  
  status(detail = "Starting Simulation", value = 0.3, session = session)
  sim <- list()
  
  for(i in samp){
    status(detail = sprintf("Running Simulation for sample %s of %s", which(i == samp),
                            length(samp)),
           value = which(i==samp)/length(samp), 
           session = NULL)
    
    sim[[paste(i)]] <- simulate_dataset(data = data, stats = stats, n_group = n_group,
                                        n_prots = n_prots, n_sim = num_simulation, 
                                        n_sample = i, exp_fc = fc,
                                        sel_sim_proteins = sel_simulated_proteins,
                                        prot_prop = prot_proportion,
                                        max_prots = prot_number, 
                                        list_diff_prots = diff_prots,
                                        sim_validation = as.logical(sim_valid),
                                        n_sample_val = valid_samples_per_grp)
  }
  
  status(detail = "Simulation Complete", value = 0.9, session = session)
  return(sim)
}


impute_random <- function(df){
  #get number of missing values of protein abundance
  n_missing <- df[is.na(abundance), .N]
  #imputed the missing values randomly from the available values
  imputed <- df[!is.na(abundance), sample(abundance, n_missing, replace = T)]
  #insert the imputed values in the data.table 
  df[is.na(abundance), abundance := imputed]
  return(df)
}

sample_simulation <- function(c_stat, n_sample, n_group){
  #simulate normally distributed data given the mean and standard deviation
  #for each protein and condition in the data
  sim <- c_stat[, rnorm(n_sample, mean, sigma), .(protein, condition)]
  sim[, ID := sequence(.N), by = .(protein,condition)]
  sim <- dcast(sim,ID+condition~protein, value.var = 'V1')
  sim[, ':='(condition = as.factor(condition), ID = NULL)]
  return(sim)
}


simulate_dataset <- function(data, stats, n_group, group, n_prots, n_sim, n_sample,
                             exp_fc = 'data', sel_sim_proteins = 'proportion',
                             prot_prop = 1, max_prots = 1000, list_diff_prots = NULL, 
                             sim_validation = FALSE, create_val = FALSE,
                             n_sample_val =NULL){
  
  if(is.null(data) || is.null(stats)){
    stop(Sys.time()," : No Data and Summary Statistics Provided")
  }
  
  #identify the number and proportion of the protiens to be simulated
  prot_num <- ifelse(sel_sim_proteins == "proportion", round(n_prots * prot_prop),
                     max_prots)
  #subset the required number of proteins
  if(!is.null(list_diff_prots)){
    stats <- stats[protein %in% list_diff_prots]
  }
  
  stats <- stats[,.SD[order(sampleMean, decreasing = T)][1:prot_num]]
  stats <- stats[!is.na(protein)]
  stats[, sampleMean := NULL]
  #convert data to long form
  stats <- melt(stats, id.vars = 'protein')
  #split the condition and statistics 
  stats[, c("stat", "condition") := tstrsplit(variable, "_.", fixed=TRUE)]
  stats <- dcast(stats, protein+condition~stat, value.var = 'value')
  
  #identify the fold change and apply it to the mean and standard deviations
  if(exp_fc != 'data'){
    exp_fc[, condition := gsub(' ','.', condition)] #check for space in conditions
    stats <- merge(stats, exp_fc, by = 'condition') #merge with mean table
    stats[, mean := fc_val * mean] # apply fold change values
    stats[, fc_val := NULL] #remove the unwanted column
  }
  
  #simulated validation data if required
  if(sim_validation  && create_val){
    message(Sys.time()," Validation Set Simulated")
    validation_data <- sample_simulation(c_stat = stats, n_sample = n_sample_val,
                                         n_group = n_group)
  }else{
    df_list <- split(data, by = 'protein')
    message(Sys.time()," Validation Set Checked for Missing Values")
    #else check if the data has any missing value to be imputed  for validation
    dest_frame <- sprintf("Sample_%s_validation", n_sample)
    validation_data <- rbindlist(lapply(df_list, impute_random))
    validation_data <- dcast(validation_data, bioreplicate+condition~protein,
                             value.var = 'abundance')
    validation_data[, bioreplicate := NULL]
    assign(dest_frame, h2o::as.h2o(validation_data))
  }
  
  simulated <- list()
  for(i in seq_len(n_sim)){
    status(detail = sprintf("Simulating ",i," of ", n_sim),
           )
    sim_data <- sample_simulation(c_stat = stats, 
                                  n_sample = n_sample,
                                  n_group = n_group)
    dest_frame_tr <- sprintf("Sample_%s_Simulation_%s", n_sample, i)
    simulated[[dest_frame_tr]] <- h2o::as.h2o(sim_data, dest_frame_tr)
  }
  message(Sys.time()," Simulation Completed")
  
  return(list('simulated' = simulated, 'validation' = get(dest_frame)))
}

#### Classification #####

sample_size_classification <- function(n_samp, sim_data, classifier, par = T, session = NULL){
  samp <- unlist(strsplit(n_samp,','))
  df <- data.table(Parameter = c("Sample Size", "Pred Accuracy", "Feature Importance"),
                   Value = c(NA, NA, NA))
  resi <- f_imp <- pred_acc <- list()
  if(classifier == 'nnet')
    par <- F
  
  for(i in seq_along(samp)){
    val <- i/length(samp)
    status(detail = sprintf("Classifying Sample Size %s of %s", i, length(samp)),
           session = session, value = val)
    res <- MSstatsSampleSize::designSampleSizeClassification(
      simulations = sim_data[[samp[i]]], classifier = classifier, parallel = par 
    )
    resi[[as.character(samp[i])]] <- res
    f_imp[[as.character(samp[i])]] <- res$feature_importance
    pred_acc[[as.character(samp[i])]] <- res$predictive_accuracy
  }
  
  return(list('res' = resi, 'samp' = as.numeric(samp), 'pred_acc' = pred_acc,
              'f_imp' = f_imp))
  }


####### H2o Application for classification #######

ss_classify_h2o <- function(n_samp, sim_data, classifier, stopping_metric = "AUTO",
                            seed = -1, nfolds = 0, fold_assignment = "AUTO", iters = 200,
                            alpha = 0, family, solver, link, min_sdev, laplace, eps,
                            session = NULL){
  samp <- unlist(strsplit(n_samp,','))
  max_val <- 0
  iter <- 0
  
  modelz <- list()
  for(i in samp){
    
    train_list <- sim_data[[i]]$simulated
    
    for(index in seq_along(train_list)){
      if(max_val == 0){
        max_val <- length(train_list) * length(samp)
      }
      iter = iter + 1/max_val
      
      status(detail = sprintf("Classifying Sample Size %s of %s, Simulation %s of %s",
                              which(samp == i), length(samp), index, length(train_list)),
             session = session, value = iter)
      
      
      if(classifier == "rf"){
        model <- h2o::h2o.randomForest(y = 1, training_frame = train_list[[index]],
                                       stopping_metric = stopping_metric, seed = seed, 
                                       balance_classes = FALSE, nfolds = nfolds,
                                       fold_assignment = fold_assignment)
        var_imp <- h2o::h2o.varimp(model)
        sel_imp <- var_imp$variable[1:10]
        model <- h2o::h2o.randomForest(y = 1, 
                                       training_frame = train_list[[index]][,c('condition',
                                                                               sel_imp)],
                                       validation_frame = sim_data[[i]]$validation,
                                       stopping_metric = stopping_metric, seed = seed, 
                                       balance_classes = FALSE, nfolds = nfolds,
                                       fold_assignment = fold_assignment)
        
      } else if (classifier == "nnet"){
        l1 = 0
        l2 = 0
        rho = 0.99
        epochs = 10
        hidden = c(250,250)
        activation = "Rectifier"
        model <- h2o::h2o.deeplearning(y = 1, 
                                       training_frame = train_list[[index]],
                                       l1= l1, l2=l2,
                                       activation = activation,
                                       hidden = hidden,
                                       epochs = epochs)
        var_imp <- h2o::h2o.varimp(model)
        sel_imp <- found_imp$variable[1:10]
        model <- h2o::h2o.deeplearning(y = 1,
                                       training_frame = train_list[[index]][,c('condition',
                                                                               sel_imp)],
                                       validation_frame = sim_data[[i]]$validation,
                                       l1 = l1, l2 = l2,
                                       activation = activation, hidden = hidden,
                                       epochs = epochs)
        
      } else if (classifier == "svmLinear"){
        model <- h2o::h2o.psvm(y = 1, training_frame = train_list[[index]], 
                               max_iterations = iters,
                               seed = seed,
                               validation_frame = sim_data[[i]]$validation,
                               disable_training_metrics = F)
        var_imp <- h2o::h2o.varimp(model)
        
      } else if (classifier == "logreg"){
        model <- h2o::h2o.glm(y = 1, training_frame = train_list[[index]],
                              seed = seed, family = family, alpha = alpha, 
                              nfolds = nfolds, solver = solver,
                              link = link)
        
        var_imp <- h2o.varimp(model)
        sel_imp <- var_imp$variable[1:10]
        model <- h2o::h2o.glm(y = 1, 
                              training_frame = train_list[[index]][,c('condition',
                                                                      sel_imp)],
                              seed = seed, family = family, alpha = alpha,
                              nfolds = nfolds, solver = solver, link = link,
                              validation_frame = sim_data[[i]]$validation)
        
      } else if (classifier == "naive_bayes"){
        model <- h2o::h2o.naiveBayes(y = 1, training_frame = train_list[[index]],
                                     ignore_const_cols = TRUE,
                                     nfolds = nfolds,
                                     fold_assignment = fold_assignment,
                                     seed = seed, laplace = laplace,
                                     min_sdev = min_sdev,
                                     eps_sdev = eps, 
                                     validation_frame = sim_data[[i]]$validation)
        var_imp <- h2o::h2o.varimp(model)
      } else{
        stop("Not defined")
      }
      
      pred <- predict(model, sim_data[[i]]$validation[,-1])
      pred_vals <- as.vector(apply(pred[-1], 1, which.max))
      pred_vals <- colnames(pred)[-1][pred_vals]
      cm <- table(pred_vals, as.vector(sim_data[[i]]$validation[,1]))
      
      acc <- sum(diag(cm))/sum(cm)
      name_val <- names(train_list)[index]
      modelz[[name_val]] <- list('model'= model, 'accuracy' = acc,
                                 'var_imp' = var_imp)
    }
  }
  return(list('models' = modelz))
}

h2o_config <- function(){
  config <- list()
  config$threads <- as.numeric(Sys.getenv('nthreads'))
  config$max_mem <- NULL
  mem <- Sys.getenv("max_mem")
  if(mem!= '')
    config$max_mem <- mem
  config$log_dir <- Sys.getenv("log_dir")
  config$log_level <- Sys.getenv("log_level")
  
  config$threads <- ifelse(is.na(config$threads), -1, config$threads)
  
  config$log_dir <- ifelse(config$log_dir == "", getwd(), config$log_dir)
  config$log_level <- ifelse(config$log_level == "", "INFO", config$log_level)
  
  return(config) 
}

plot_acc <- function(data, use_h2o, alg = NA){
  if(use_h2o){
    shiny::validate(shiny::need(data$models, "No Models Run Yet"))
    #loop through the object returned by classification to extract accuracy
    model_data <- data$models
    df <- rbindlist(lapply(names(model_data), function(x){
      z <- model_data[[x]]
      strs <- unlist(strsplit(x,' '))
      
      cm <- z$model@model$validation_metrics@metrics$cm$table
      cm <- cm[1:(dim(cm)[1]-1),1:(dim(cm)[2] -2)]
      cm <- as.matrix(sapply(cm, as.numeric))
      acc <- sum(diag(cm))/sum(cm)
      
      data.table(sim  = as.numeric(gsub("[[:alpha:]]",'',strs[2])),
                 sample = as.factor(gsub("[[:alpha:]]",'',strs[1])),
                 mean_acc = acc)
    }))
  }else{
    shiny::validate(shiny::need(data$samp, "No Trained Models Found"))
    df <- rbindlist(data$pred_acc)
    names(df) <- c("sample","mean_acc")
  }
  
  df[, acc := mean(mean_acc), sample]
  setorder(df, -sample)
  # following logic is flawed only workds if the data.table is arrange by sample
  # size in increasing order
  #######
  mean_PA <- df$acc
  sample_size <- df$sample
  dydx <- -diff(mean_PA)/-diff(as.numeric(as.character(sample_size)))
  
  if(any(dydx >= 0.0001)){
    inter_dydx <- dydx[which(dydx >= 0.0001)]
    min_dydx <- inter_dydx[which.min(inter_dydx)]
    optimal_index <- which(dydx == min_dydx)
    optimal_sample_size_per_group <- sample_size[optimal_index]
  } else{
    optimal_sample_size_per_group <- sample_size[1]
  }
  
  y_lim <- c(df[,min(acc, na.rm = T)]-0.1, 1)
  df[sample == optimal_sample_size_per_group, fill_col := 'red']
  ######
  
  p <- ggplot(data = df, aes(x = reorder(sample)))+
    geom_boxplot(aes(y = mean_acc, group = sample, fill = fill_col), alpha = 0.5)+
    scale_fill_identity()+
    # geom_vline(xintercept = optimal_sample_size_per_group, color = 'red', size= 0.75)+
    geom_point(aes(y = acc))+
    geom_line(aes(y = acc, group = 1), size = 0.75, color = "blue")+
    labs(x = "Simulated Sample Size", y = "Predictive Accuracy",
         title = sprintf("Classifier %s", alg),
         subtitle = sprintf("Optimum accuracy achieved when sample size is : %s",
                            optimal_sample_size_per_group))+
    ylim(y_lim)+
    theme_MSstats()+
    theme(plot.subtitle = element_text(face = 'italic', color = 'red'))
  
  return(p)
}

plot_var_imp <- function(data, sample = 'all', sim = NA, alg = NA, use_h2o, prots = 10){
  if(use_h2o){
    if(prots == 'all'){
      prots <- nrow(data$models[[1]]$var_imp)
    }
    
    shiny::validate(shiny::need(length(data) != 0, "No Trained Models Found"))
    if(sample == 'all'){
      sample <- names(data$models)
    }else{
      sample <- gsub('Sample','',sample)
      sample <- sprintf("Sample_%s_Simulation_%s", sample, 1:as.numeric(sim))
    }
    df <- rbindlist(lapply(sample, function(x){
      dt <- as.data.table(data$models[[x]]$var_imp)
      dt$name <- x
      dt
    }))
    
    df[, c('noreq', 'sample_size', 'noreq2', 'simulation') := tstrsplit(name, "_", fixed = T)]
    df[,':='(name = NULL, noreq = NULL, noreq2 = NULL)]
    df <- split(df, by = 'sample_size')
    
    dt <- rbindlist(lapply(df, function(x){
      samp <- unique(x$sample_size)
      d <- dcast(x, variable~sample_size+simulation,
                 value.var = 'relative_importance')
      d <- cbind(d[,1], rowSums(d[,-1], na.rm = T))
      setorder(d,-V2)
      d[, ':='(V2 = (V2-min(V2))/(max(V2) - min(V2)),
               sample_size = samp)]
      d[1:prots]
    }))
    setnames(dt, 'V2', 'relative_importance')
  }else{
    if(sample == 'all'){
      sample <- as.character(data$samp)
    }else{
      sample <- gsub('Sample','',sample)
    }
    
    df <- rbindlist(lapply(sample, function(x){
      d <- data$f_imp[[x]]
      d[, sample_size := paste("SampleSize",x)]
      setnames(d, c('protein.rn', 'importance'),
               c('variable', 'relative_importance'),
               skip_absent = T)
    }))
    
    if(prots == 'all'){
      prots <- max(df[,.N,sample_size][,unique(N)], na.rm = T)
    }
    setorder(df, -relative_importance)
    dt <- df[,head(.SD, prots),,by = c("sample_size")]
  }
  
  dt <- split(dt, by = 'sample_size')
  
  g <- lapply(dt,  function(x){
    x$variable <- reorder(x$variable, x$relative_importance)
    x[, sample_size := gsub('SampleSize','', sample_size)]
    ggplot(data = x, aes(variable, relative_importance))+
      geom_col()+
      labs(x = "Protein", y = "Relative Importance", 
           title = paste('Sample Size', unique(x$sample_size)))+
      scale_x_discrete(breaks = x$variable,
                       labels = gsub("_.*",'',as.character(x$variable)))+
      theme_MSstats()+
      coord_flip()
  })
  names(g) <- names(dt)
  return(g)
}

#### WRAPPER FOR CLASSIFICATION ######

run_classification <- function(sim, inputs, use_h2o, seed, session = session){
  if(seed != -1)
    set.seed(seed)
  
  if(use_h2o){
    classification <- ss_classify_h2o(n_samp = inputs$n_samp_grp, sim_data = sim,
                                      classifier = inputs$classifier,
                                      stopping_metric = inputs$stop_metric,
                                      nfolds = inputs$nfolds,
                                      fold_assignment = inputs$f_assignment, iters = inputs$iters,
                                      family = inputs$family, solver = inputs$solver,
                                      link = inputs$link, min_sdev = inputs$min_sdev,
                                      laplace = inputs$laplace, eps = inputs$eps_sdev,
                                      seed = -1, session = session)
  }else{
    
    classification <- sample_size_classification(n_samp = inputs$n_samp_grp,
                                                 sim_data = sim,
                                                 classifier = inputs$classifier,
                                                 session = session)
  }
  return(classification)
}

sample_size_classification <- function(n_samp, sim_data, classifier, k = 10,
                                       family = 'binomial', session = NULL){
  samp <- unlist(strsplit(n_samp,','))
  pred_acc <- list()
  f_imp <- list()
  models <- list()
  max_val <- 0
  iter <- 0
  
  for(i in seq_along(samp)){
    res <- list()
    imp <- list()
    model <- list()
    list_train <- sim_data[[i]]$simulated
    
    fam_check <- unique(as.vector(sim_data[[i]]$validation[, 1]))
    
    if(length(fam_check) > 2){
      family <- 'multinomial'
    }
    
    tryCatch({
      for(j in seq_along(list_train)){
        if(max_val == 0){
          max_val <- length(list_train) * length(samp)
        }
        iter = iter + 1/max_val
        
        status(detail = sprintf("Classifying Sample Size %s of %s, Simulation %s of %s",
                                i, length(samp), j, length(list_train)),
               session = session, value = iter)
        
        res[[paste0('Sim',j)]] <- classify(df = as.data.table(list_train[[j]]), 
                                           val = as.data.table(sim_data[[i]]$validation), 
                                           alg = classifier,
                                           family = family, k = k)
      }
    }, error = function(e){
      print(e)
    })
    
    for(j in seq_along(res)){
      acc <- data.frame('Sample' = samp[i],'accuracy' = res[[j]]$accuracy)
      pred_acc <- append(list(acc), pred_acc)
      
      imp[[j]] <- data.table('Simulation' = j, res[[j]]$f_imp[1:k])
      model[[j]] <- res[[j]]$model
    }
    
    imp <- do.call('rbind', imp)
    imp <- dcast(imp, rn~Simulation, value.var = 'Overall')
    imp <- cbind('protein' = imp[,1],
                 'importance' = rowSums(imp[,-1], na.rm = T))
    imp[, importance := (importance-min(importance))/(max(importance) - min(importance))]
    
    models[[as.character(samp[i])]] <- model
    f_imp[[as.character(samp[i])]] <- imp
  }
  return(list('res' = models, 'samp' = as.numeric(samp), 'pred_acc' = pred_acc,
              'f_imp' = f_imp))
}



classify <- function(df, val, alg, family, k){
  if(alg == 'logreg'){
    alg = 'glm'
  }
  
  if(alg == 'glm'){
    if(family == 'multinomial'){
      model <- nnet::multinom(condition~., data =df , maxit, MaxNWts = 84581)
      f_imp <- caret::varImp(model, scale = T)
      sel_imp <- rownames(f_imp)[1:k]
      sel_imp <- gsub('`','',sel_imp)
      if(!all(sel_imp %in% names(df))){
        sel_imp <- gsub('`','',sel_imp)
      }
      model <- nnet::multinom(condition~., data = df[,, c('condition', sel_imp)],
                              maxit=1000,MaxNWts=84581)
    } else {
      model <- caret::train(make.names(condition)~.,data = df,
                            method = alg, family = family,
                            trControl = caret::trainControl(method = "none",
                                                            classProbs = TRUE))
      f_imp <- caret::varImp(model, scale = TRUE)
      i_ff <- data.table::as.data.table(f_imp$importance, keep.rownames = T)
      setorder(i_ff, -Overall)
      sel_imp <- i_ff[1:k, rn]
      if(!all(sel_imp %in% names(df))){
        sel_imp <- gsub('`','',sel_imp)
      }
      model <- caret::train(make.names(condition)~., 
                            data = df[,, c('condition', sel_imp)], 
                            method = alg,
                            trControl = caret::trainControl(method = "none",
                                                            classProbs = TRUE))
    }
  } else {
    model <- caret::train(make.names(condition)~., data = df,
                          method = alg, 
                          trControl = caret::trainControl(method = "none", 
                                                          classProbs = TRUE)) 
    
    f_imp <- caret::varImp(model, scale = TRUE)
    i_ff <- data.table::as.data.table(f_imp$importance, keep.rownames = T)
    setorder(i_ff, -Overall)
    sel_imp <- i_ff[1:k, rn]
    
    if(!all(sel_imp %in% names(df))){
      sel_imp <- gsub('`','',sel_imp)
    }
    model <- caret::train(make.names(condition)~.,
                          data = df[,, c('condition', sel_imp)],
                          method = alg, 
                          trControl = caret::trainControl(method = "none", 
                                                          classProbs = TRUE)) 
    
  }
  
  pred <- predict(model, val[,-1])
  cm <- table(pred, val$condition)
  
  acc <- sum(diag(cm))/sum(cm)
  r_list <- list('accuracy' = acc, 'model' = model, 'f_imp' = i_ff)
  return(r_list)
}




