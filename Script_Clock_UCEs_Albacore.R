#############################################################################################################################################################
##################################################################                                ###########################################################
##################################################################        MammalMethyl Clock      ###########################################################
##################################################################                                ###########################################################
#############################################################################################################################################################
library(data.table)
library(stringr)
library(MammalMethylClock)
library(dplyr)
library(tidyverse)
library(ggplot2)

#### Path #### 
root_dir   <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/"
setwd(root_dir)
meth_file  <- file.path(root_dir, "Metadata_age_tuna.csv")
meta_file  <- file.path(root_dir, "Jesstimation_ALB_papier.csv")

meth <- fread(meth_file) |> as.data.frame()
stopifnot("cpg_id" %in% names(meth))
meth <- meth[,-1]
rownames(meth) <- meth$cpg_id
meth <- meth[,-1]

meta <- fread(meta_file) |> as.data.frame()
stopifnot(all(c("sample_id","age") %in% names(meta)))
meta$age <- as.numeric(meta$age)

meth_matrix <- t(meth)
age_tuna_final <- meta 

#### Clock #### 

OUTVAR <- "age"
PREDVAR <- "DNAmAgebasedOnAll"
RESVAR <- "AgeAccelbasedOnAll"
RESinOtherVAR <- "AgeAccelinOtherbasedOnAll"
out.csv <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/Tuna_Clock_Final_EpigeneticLLin3Age_Log_transform.csv"
output.csv <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/Tuna_Clock_Final_EpigeneticLLin3Age_Log_transform_PredictedValues.csv"
out.png <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/Tuna_Clock_Final_EpigeneticLLin3Age_Log_transform.tiff"
out.png.title <- "Tuna_Clock_basedOnAll_EpigeneticLLin3Age"
fun_VAR1 <- NULL
fun_VAR2 <- NULL
ALPHA <- 0.5
NFOLD <- 10
ys.output <- saveBuildClock(
  xs.train = meth_matrix, ys.train = age_tuna_final, xs.other = NULL, ys.other = NULL, OUTVAR,
  out.csv, output.csv, out.png, out.png.title,
  PREDVAR, RESVAR, RESinOtherVAR, ALPHA, NFOLD,
  fun_trans=fun_log.trans,fun_inv=fun_log.inv
  # ,  fun_trans = fun_llin3.trans, fun_inv = fun_llin3.inv,
  # fun_VAR1 = fun_VAR1, fun_VAR2 = fun_VAR2
)


#### LOOCV #### 

OUTVAR <- "age"
COLVAR <- "Tissue"
PREDVAR <- "DNAmAgeLOO"
RESVAR <- "AgeAccelLOO"
out.rdata <- NULL
output.csv <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/Tuna_LOO_Final_EpigeneticLLin3Age_Log_transform.csv"
out.png <- "D:/Analyse_Thomas/MammalMethyl/Tuna_UCEs/Tuna_LOO_Final_EpigeneticLLin3Age_Log_transform.png"
out.png.title <- "Tuna_LOO_Final_EpigeneticLLin3Age"
fun_VAR1 <- NULL
fun_VAR2 <- NULL
ALPHA <- 0.5
NFOLD <- 10
ys.output <- saveLOOEstimation(
  xs = meth_matrix, ys = age_tuna_final, OUTVAR,
  out.rdata, output.csv, out.png, out.png.title,
  PREDVAR, RESVAR, ALPHA, NFOLD,
  fun_trans=fun_log.trans,fun_inv=fun_log.inv
  # , fun_trans = fun_llin3.trans, fun_inv = fun_llin3.inv,
  # fun_VAR1 = fun_VAR1, fun_VAR2 = fun_VAR2, COLVAR = COLVAR
)


















