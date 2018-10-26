###---------------------------------------------------------------------------------------------------------------------------
# 
#     LUNG SURVIVAL CLASSIFICATION (DEEP LEARNING Vs. ENGINEERED FEATURES)
#     Initial edit: 05/17/17
#     Last Edit: 03/27/18
#     tcoroller@lroc.harvard.edu
#     
###---------------------------------------------------------------------------------------------------------------------------
rm(list = ls()); graphics.off() # Clean environment


###---------------------------------------------------------------------------------------------------------------------------
### Setup
###--------------------------------------------------------------------------------------------------------------------------

# set path
setwd('Z:/USERS/Ahmed/LUNGS/Plos Medicine - Submission/') # Root files

# Load radR package
# devtools::install_github("CompRad/radR", ref = '3.0', auth_token = "d0435d8038e71e6dee305f843bd3e88e92d560b8") # For first use only

require(dplyr)
require(tidyr)
require(reshape2)
require(pROC)
require(pheatmap)
require(survcomp)
require(ggplot2)
require(RColorBrewer)
require(radR)
require(irr)


# Main function made for the analysis
source("./2_R_analysis/Function_models.R")


###---------------------------------------------------------------------------------------------------------------------------
### Curation
###--------------------------------------------------------------------------------------------------------------------------

# Netherlands
maastro <- read.csv(file = "./0_data/Maastro.csv")
radbound <- read.csv(file = "./0_data/Radboud.csv")
mumc <- read.csv(file = "./0_data/MUMC.csv")

# Moffit
moffitt <- read.csv(file = "./0_data/Moffitt.csv")
m_spore <- read.csv(file = "./0_data/M-SPORE.csv")

# Harvard
harvardRT <- read.csv(file = "./0_data/HarvardRT.csv")


# 2029 features
engineeredColumns <- colnames(harvardRT)[grep(pattern = "original|wavelet|log|square|expo", colnames(harvardRT))]
engineeredColumns <- engineeredColumns[-grep(pattern = "histology|logit|stage", engineeredColumns)]


# Excluding patients with NAs
maastro <- maastro[-grep(pattern = 'LUNG1-286', x = maastro$id),]
harvardRT <- harvardRT[!harvardRT$id == 151,]
harvardRT <- harvardRT[!harvardRT$id == 241,]

# Convert back to numeric - Issue with the patients with NAs that converts the features as character.
maastro[,engineeredColumns] <- as.data.frame(sapply(X = engineeredColumns, function(X) as.numeric(as.character(maastro[,X]))))
harvardRT[,engineeredColumns] <- as.data.frame(sapply(X = engineeredColumns, function(X) as.numeric(as.character(harvardRT[,X]))))



###==================================
###    ** Stability Assessment **
###==================================


# Load data
rider_test <- read.csv('./0_data/RIDER_Test.csv')
rider_retest <- read.csv('./0_data/RIDER_Retest.csv')


## 2 - compute ICC
ICC <- sapply(X=intersect(engineeredColumns, colnames(rider_test)), function(X){
  abs(irr::icc(cbind(rider_test[,X], rider_retest[,X]), model="twoway", type = "a")$value)
})

ICC[is.na(ICC)] <- 0 # safety step

# Keep features with ICC above 0.8, following guidelines Parmar et al, 2015
f.robust <- engineeredColumns[ICC > 0.8]

# remove anything with NAs, safety step
f.robust <- setdiff(f.robust, colnames(rider_retest)[ apply(rider_retest, 2, anyNA) ])




###===================================
###    ** Prep for classification **
###===================================


f.engineered <-  intersect(intersect(intersect(colnames(moffitt), f.robust), colnames(mumc)), colnames(harvardRT)) # Safety check
f.clinical <- c('age','stage','histology_grouped', "gender")
f.deep <- "logit_1"
Class <- "surv2yr"

#----------------------------------------
# ** Curation 
maastro$Cohort <- "maastro"
radbound$Cohort <- "radbound"
mumc$Cohort <- "mumc"
harvardRT$Cohort <- "harvard-rt"
moffitt$Cohort <- "moffitt"
m_spore$Cohort <- "m-spore"

# Combine dataset
f.tokeep <- c(Class, f.clinical, f.engineered, f.deep, "Cohort", "deadstat")
df <- dplyr::bind_rows(maastro[,f.tokeep], 
                       radbound[,f.tokeep], 
                       mumc[,f.tokeep], 
                       moffitt[,f.tokeep], 
                       m_spore[,f.tokeep], 
                       harvardRT[,f.tokeep], .id = "id")
        
        

# Factorise Class for classifiers - need levels
df[,Class] <- paste("class", df[,Class], sep = ".") %>% as.factor()

# Remove NaN
f.engineered <- f.engineered[colSums(as.data.frame(apply(df[,f.engineered], 2, is.nan)))==0]
df <- df[,f.tokeep <- c(Class, f.clinical, f.engineered, f.deep, "Cohort", "deadstat")]




#---------------------------------------------------
# ** RUN ANALYSIS **
#---------------------------------------------------


f.volume <- c("original_shape_Volume", "original_shape_Maximum2DDiameterSlice")


# RT
list.RT <- Function_models(Class = Class, df = df, select = "radbound", train =  "harvard-rt", test = "maastro",
                           f.engineered = f.engineered, f.volume = f.volume, f.clinical = c('age','stage','histology_grouped'))

# SURGERY
list.Surgery <- Function_models(Class = Class, df = df, select = "mumc", train = "moffit", test = "m-spore",
                                f.engineered = f.engineered, f.volume = f.volume, f.clinical = c('gender','stage','histology_grouped'))


# CURATION
results.RT <- list.RT$d.results
results.Surgery <- list.Surgery$d.results


cat("\n*************************\n")
print(results.RT)
print(results.Surgery)
cat("\n*************************\n")


results.Surgery$Type <- "Surgery"
results.RT$Type <- "RT"


d.plot <- dplyr::inner_join(reshape2::melt(rbind(results.RT[1,], results.Surgery[1,]))  %>%
                            dplyr::rename(AUC = value) %>%
                            dplyr::select(variable, AUC, Type),
                            reshape2::melt(rbind(results.RT[2,], results.Surgery[2,]))  %>%
                            dplyr::rename(p.value = value) %>%
                            dplyr::select(variable, p.value, Type)) %>%
                            select(variable, AUC, p.value, Type)


d.plot <- dplyr::inner_join(reshape2::melt(rbind(results.RT[3,], results.Surgery[3,]))  %>%
                            dplyr::rename(permutation = value) %>%
                            dplyr::select(variable, permutation, Type), 
                            d.plot)


d.plot$variable <- factor(x = as.character(d.plot$variable), 
                          levels = c('Deep', 'Engineered', 'Deep_Engineered', 'Clinical',
                                     'Deep_Clinical', 'Volume', 'Diameter', 'All', 'Engineered_Clinical'))


#---------------------------------------------------
# ** PLOT **
#---------------------------------------------------


plot_auc <- ggplot(data = d.plot, mapping = aes(x = variable, y = AUC)) +
    
    geom_bar(stat="identity", aes(fill=variable)) + 
    
    theme_bw() +
    
    theme(legend.position = "none",
          axis.text.x=element_text(angle = 45, hjust = 1)) +
    
    labs(x = "Models", y = "AUC") +
    
    coord_cartesian(ylim = c(0.5, 0.75)) +
    
    facet_grid(facets = Type ~ . ) +
  
    geom_text(aes(label = format(AUC, digits = 2, scientific = F)), size = 4, vjust = 1.2, fontface = "bold") +
    
    geom_text(data = filter(d.plot, permutation < 0.05), label = ">", size = 7, vjust = -0.5) +
    
    geom_text(data = filter(d.plot, variable == "Deep"), label = "ref", size = 5, vjust = -0.5) +
    
    
    geom_text(data = filter(d.plot, p.value < 0.05), label = "*", size = 10, y = 0.53, color = 'gray10') +
    
    geom_text(data = filter(d.plot, p.value < 0.05), aes(label = format(p.value, digits = 3, scientific = T)), size = 4, y = 0.51)
    



#---------------------------------------------------
# ** Meta p-values **
#---------------------------------------------------

models <- c("o.eng", "o.deepEng", "o.DeepClin", "o.EngClin", "o.All", "o.clin", "o.volume", "o.diameter")

meta.pval <- sapply(X =models , function(X){
  
  cindex.comp.meta(list.cindex1 = list("deep.RT" = list.RT$o.deep, "Deep.Sx" = list.Surgery$o.deep), 
                   list.cindex2 = list("X.RT" = list.RT[[X]], "X.Sx" = list.Surgery[[X]]), hetero = T)$p.value
  
})

meta.pval <- data.frame(models = paste("Deep", models, sep = " + "),
                        meta.pvals = as.numeric(meta.pval))


###==================================
###      ** EXPORT ALL  **
###==================================


RESULTS <- list(RT = list.RT,
                Surgery = list.Surgery, 
                TABLE_AUC = d.plot,
                TABLE_meta_p = meta.pval, 
                dataset = df,
                f.robust = f.robust,
                f.engineered = f.engineered,
                plot_auc = plot_auc)

filename <- paste('OUT_results_', substr(Sys.time(), start = 1, stop = 10), '.rda', sep = "")
save(RESULTS, file = file.path("./2_R_analysis/", filename))

