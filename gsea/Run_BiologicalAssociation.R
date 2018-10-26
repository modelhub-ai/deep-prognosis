############################################################################
###
### Associate deep learned vectors with biological pathways.
###
### Patrick Grossmann, June 12, 2017, patrickb_grossmann@dfci.harvard.edu
###
############################################################################

###
### This script will calculate enrichment statistics for two datasets
###

# ----------------------------------------------------------------------#
# ------------------- Read in deep learning vectors ------------------- #
# ----------------------------------------------------------------------#

getVecAsTab <- function(fn) {
  tab <- read.csv(fn)
  vec <- tab[,"logit_1"]
  names(vec) <- tab[,"patient"]
  data.frame(imaging.DL=vec)
}

DL.moffitt <- getVecAsTab("moffitt_train.csv")
DL.lung3 <- getVecAsTab("lung3_test.csv")

# ----------------------------------------------------------------#
# ------------------- Read in biological data ------------------- #
# ----------------------------------------------------------------#

load("moffitt.curation.rda")
expr.moffitt <- moffitt.uncombined$expression
load("maastro.curation.rda")
expr.lung3 <- maastro.uncombined$expression

# ---------------------------------------------------#
# ------------------- Rank genes ------------------- #
# ---------------------------------------------------#
require(RadioGx)

rank.moffitt <- makeRanking2(geno=expr.moffitt, pheno=DL.moffitt, outdir="ranking", name="Moffitt", 
                             pref.pdata="imaging\\.", pref.exprs="geneid\\.", cores=1)
rank.lung3 <- makeRanking2(geno=expr.lung3, pheno=DL.lung3, outdir="ranking", name="Lung3", 
                           pref.pdata="imaging\\.", pref.exprs="geneid\\.", cores=1)

# --------------------------------------------------------------#
# ------------------- Associate to pathways ------------------- #
# --------------------------------------------------------------#

genesets <- "~/research/genesets/c2.cp.reactome.v6.0.entrez.gmt"
algorithm <- "~/research/algorithm/gsea-3.0.jar"

gsea.Moffitt <- performGSEA(dsName="Moffitt", ranks=rank.moffitt$files[["imaging.DL"]],
                            min.size=5, max.size=500, gmt.path=genesets, 
                            out.dir="GSEA/Moffitt", exec.path=algorithm, lrdataDir="rdata", 
                            lintermediateDir="rdata/intermediate", 
                            seed=987654321, nPerm=1000, cores=1, verbose=F,
                            reuseSaved=F) 

gsea.Lung3 <- performGSEA(dsName="Lung3", ranks=rank.lung3$files[["imaging.DL"]],
                          min.size=5, max.size=500, gmt.path=genesets, 
                          out.dir="GSEA/Lung3", exec.path=algorithm, lrdataDir="rdata", 
                          lintermediateDir="rdata/intermediate", 
                          seed=987654321, nPerm=1000, cores=1, verbose=F,
                          reuseSaved=F) 

## save results
save(gsea.Moffitt, gsea.Lung3, file="gsea_results.rda")


