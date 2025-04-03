############################################ 1. Initialisation Section ############################################
rm(list=ls())
#---------------------------------- 1.1 Set up package libraries ----------------------------------
# List of packages we need
list.of.packages = c("tibble", "smoof", "flacco", "lhs", "RANN", "devtools", "numDeriv", "e1071", "expm", "mda",
                     "mlbench", "shape", "testthat", "openxlsx", "nsga2R", "R.matlab", "gtools", "dplyr", "corrplot",
                     "stringi", "stringr")

# Install packages if they are not already installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
for (i in 1:length(list.of.packages)){
  library(list.of.packages[i], character.only=TRUE, quietly=TRUE)
}

featurelist <- listAvailableFeatureSets(allow.additional_costs = FALSE)[-c(1,2,3,12,13)]

# Set working directory
array.idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (array.idx>5400){
  indir <- '~/punim0320/CORNN_data/'
  array.idx <- array.idx - 5400
} else {
  indir <- '~/punim0320/BBOB_data/'
}
outdir <- '~/punim0320/FLACCO_data/'

filelist <- list.files(indir, pattern='F')

currfile <- filelist[array.idx]

if(file.exists(paste0(outdir,"ELA_",currfile))){
  stop()
}

idx <- str_locate(currfile, ".csv")[1]-1
R <- as.numeric(str_sub(currfile, idx, idx))
if (is.na(R)){
  print(paste0(currfile,' failed due to incorrect repeat number'))
  stop()
}

print(paste0(currfile,' started!'))
Y <- read.csv(paste0(indir, currfile))
D <- nrow(Y)/100
X <- read.csv(paste0('~/punim0320/CORNN_data/X_D', as.character(D), '_S100_R', as.character(R), '.csv' ))

features <- NULL

for (jj in 2:ncol(Y)){
  print(paste0('Processing column: ',as.character(jj)))
  feat.object = createFeatureObject(X = X, y = scale(Y[,jj]))
  feat.inst <- NULL
  for (ii in 1:length(featurelist)){
    print(paste0('     -> Feature ',featurelist[ii],' started.'))
    feat <- calculateFeatureSet(feat.object, featurelist[ii])
    feat.inst <- append(feat.inst, feat)
  }
  features_append <- data.frame(t(sapply(feat.inst,c)))
  fcname <- str_remove(currfile,".csv")
  if (ncol(Y)>2){
    fcname <- paste0(fcname,"_",colnames(Y)[jj])
  }
  features <- rbind(features, cbind(fcname,features_append))
}
print(paste0(fcname,' completed'))

# Write table
write.csv(features, paste0(outdir,"ELA_",currfile)) 
