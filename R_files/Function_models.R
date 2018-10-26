


Function_models <- function(Class, df, select, train, test, classifier = "rf", nfeatures = 40,
                    f.clinical, f.engineered, f.volume, f.deep = "logit_1"){
  
  require(dplyr)
  
  number = 5
  repeats = 5
  
  # Initialize nodes
  foreach::registerDoSEQ()
  
  #---------------------------------------------------
  # ** Clinical model **
  #---------------------------------------------------
  cat("Computing Clinical...\n")
  
  # Subset
  d.train <- df %>%
    filter(grepl(train, Cohort)) %>%
    select(one_of(c(f.clinical, Class)))
  
  d.test <- df %>%
    filter(grepl(test, Cohort)) %>%
    select(one_of(c(f.clinical, Class)))

  # Run RF clinical
  set.seed(seed = 123, kind = "L'Ecuyer-CMRG")
  fit.Clin <- radR::classifier.switch(Class = Class, train = d.train, test = d.test, 
                                      classifier = classifier, number = number, repeats = repeats)
  
  #---------------------------------------------------
  # ** Engineered model **
  #---------------------------------------------------
  cat("Computing Engineered...\n")
  
  # Subset for Train / Valid
  d.select <- df %>%
    filter(grepl(select, Cohort)) %>%
    select(one_of(c(f.engineered, Class)))
  
  # mRMR
  d.select$TOKEN <- ordered(match(x = d.select[,Class], table = unique(d.select[,Class])))
  f.selection <- f.engineered
  
  # f.selection <- as.character(select.var(df = d.select, features = f.selection, corr.max = 0.95, var.number = 200))
  f.selection <- as.character(select.mrmr(train = d.select, features = f.engineered, Class = "TOKEN",
                                          method = "classic", nfeatures = nfeatures)$final.model)
  
  # Subset 
  d.train <- df %>%
    filter(grepl(train, Cohort)) %>%
    select(one_of(c(f.selection, Class)))
  
  d.test <- df %>%
    filter(grepl(test, Cohort)) %>%
    select(one_of(c(f.selection, Class)))
  
  # Initialize parallel computing
  foreach::registerDoSEQ() # Initialize nodes
  
  # Fit RF Engineered
  set.seed(seed = 123, kind = "L'Ecuyer-CMRG")
  fit.Eng <- radR::classifier.switch(Class = Class, train = d.train, test = d.test, 
                                     classifier = classifier, number = number, repeats = repeats)
  
  
  #---------------------------------------------------
  # ** Curate output **
  #---------------------------------------------------
  cat("Computing Curation...\n")
  
  scores.train <- data.frame(Clinical = fit.Clin$pred$Train$class.1,
                             Engineered = fit.Eng$pred$Train$class.1, 
                             Deep = df %>%
                               filter(grepl(train, Cohort)) %>%
                               select(one_of(f.deep))  %>%
                               dplyr::rename(Deep = logit_1),
                             True = fit.Clin$pred$Train$obs)
  
  
  scores.valid <- data.frame(Clinical = fit.Clin$pred$Valid$class.1,
                             Engineered = fit.Eng$pred$Valid$class.1, 
                             Deep = df %>%
                               filter(grepl(test, Cohort)) %>%
                               select(one_of(f.deep))  %>%
                               dplyr::rename(Deep = logit_1),
                             True = fit.Clin$pred$Valid$obs)

  
  #---------------------------------------------------
  # ** Combination **
  #---------------------------------------------------
  cat("Computing Combination...\n")
  
  
  classifier = "glm"
  # Engineered - Deep
  fit.EngDeep <- classifier.switch(Class = "True", 
                                   train = select(scores.train, Deep,  Engineered, True), 
                                   test = select(scores.valid, Deep,  Engineered, True), 
                                   classifier = classifier, number = number, repeats = repeats)
  
  # Engineered - Deep - Clinical
  fit.All <- classifier.switch(Class = "True", 
                               train = scores.train, 
                               test = scores.valid,
                               classifier = classifier, number = number, repeats = repeats)

  # Deep - Clinical
  fit.DeepClin <- classifier.switch(Class = "True", 
                                    train = select(scores.train, Deep,  Clinical, True), 
                                    test = select(scores.valid, Deep,  Clinical, True), 
                                    classifier = classifier, number = number, repeats = repeats)

  # Engineered - Clinical
  fit.EngClin <- classifier.switch(Class = "True", 
                                   train = select(scores.train, Clinical,  Engineered, True), 
                                   test = select(scores.valid, Clinical,  Engineered, True), 
                                   classifier = classifier, number = number, repeats = repeats)
 
  
  #---------------------------------------------------
  # ** Concordance Object **
  #---------------------------------------------------

  o.deep <- concordance.index(x = scores.valid$Deep, cl = scores.valid$True, method = "noether")
  o.eng <- concordance.index(x = scores.valid$Engineered, cl = scores.valid$True, method = "noether")
  o.clin <- concordance.index(x = scores.valid$Clinical, cl = scores.valid$True, method = "noether")
  o.deepEng <- concordance.index(x = fit.EngDeep$pred$Valid$class.1, cl = fit.EngDeep$pred$Valid$obs, method = "noether")
  o.All <- concordance.index(x = fit.All$pred$Valid$class.1, cl = fit.All$pred$Valid$obs, method = "noether")
  o.DeepClin <- concordance.index(x = fit.DeepClin$pred$Valid$class.1, cl = fit.DeepClin$pred$Valid$obs, method = "noether")
  o.EngClin <- concordance.index(x = fit.EngClin$pred$Valid$class.1, cl = fit.EngClin$pred$Valid$obs, method = "noether")
  
  
  d.test <- df %>% filter(grepl(test, Cohort))
  
  o.volume <- concordance.index(x = d.test[,"original_shape_Volume"], cl = d.test[,"surv2yr"], method = "noether")
  o.diameter <- concordance.index(x = d.test[,"original_shape_Maximum2DDiameterSlice"], cl = d.test[,"surv2yr"], method = "noether")
  
  
  #---------------------------------------------------
  # ** OUTPUT **
  #---------------------------------------------------
  cat("Computing Output...\n")

  # Summarize
  d.AUC <- data.frame(Deep = o.deep$c.index,
                      Engineered = o.eng$c.index,
                      Deep_Engineered = o.deepEng$c.index,
                      Deep_Clinical = o.DeepClin$c.index,
                      Engineered_Clinical = o.EngClin$c.index,
                      All = o.All$c.index,
                      Clinical = o.clin$c.index,
                      Volume = o.volume$c.index,
                      Diameter = o.diameter$c.index)
  

  d.AUC[d.AUC < 0.5] <- 1 - d.AUC[d.AUC <0.5] # Flip values
  d.AUC$Metric = "AUC"
  
  d.pvalue <- data.frame(Metric = "p.value",
                         Deep = o.deep$p.value,
                         Engineered = o.eng$p.value,
                         Deep_Engineered = o.deepEng$p.value,
                         Deep_Clinical = o.DeepClin$p.value,
                         Engineered_Clinical = o.EngClin$p.value,
                         All = o.All$p.value,
                         Clinical = o.clin$p.value,
                         Volume = o.volume$p.value,
                         Diameter = o.diameter$p.value)
  
  cat("Computing Permutations...\n")
  set.seed(1234)
  N = 1000
  sided = "lesser"
  
  d.permuts <- data.frame(Metric = "permutation",
                          Deep = 1,
                          Engineered = compare.permutation(SO1 = o.deep, SO2 = o.eng, N = N, sided = sided)$p.value,
                          Deep_Engineered = compare.permutation(SO1 = o.deep, SO2 = o.deepEng, N = N, sided = sided)$p.value,
                          Deep_Clinical = compare.permutation(SO1 = o.deep, SO2 = o.DeepClin, N = N, sided = sided)$p.value,
                          Engineered_Clinical = compare.permutation(SO1 = o.deep, SO2 = o.EngClin, N = N, sided = sided)$p.value,
                          All = compare.permutation(SO1 = o.deep, SO2 = o.All, N = N, sided = sided)$p.value,
                          Clinical = compare.permutation(SO1 = o.deep, SO2 = o.clin, N = N, sided = sided)$p.value,
                          Volume = compare.permutation(SO1 = o.deep, SO2 = o.volume, N = N, sided = sided)$p.value,
                          Diameter = compare.permutation(SO1 = o.deep, SO2 = o.diameter, N = N, sided = sided)$p.value)
  
                         
  d.results <- rbind(d.AUC, d.pvalue, d.permuts)

  d.exp <- list(d.results = d.results,
                d.permuts = d.permuts,
                fit.Clin = fit.Clin,
                fit.Eng = fit.Eng,
                fit.EngDeep = fit.EngDeep,
                fit.DeepClin = fit.DeepClin,
                fit.All = fit.All,
                fit.EngClin = fit.EngClin,
                o.deep = o.deep,
                o.eng = o.eng,
                o.deepEng = o.deepEng,
                o.DeepClin = o.DeepClin,
                o.EngClin = o.EngClin,
                o.All = o.All,
                o.clin = o.clin,
                o.volume = o.volume,
                o.diameter = o.diameter,
                f.selection = f.selection)
  
  return(d.exp)
  
  
}