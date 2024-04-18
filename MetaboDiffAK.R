library(MetaboDiff)

##First step is to download the data; raw values of metabolites as assay (matrix),
#    the treatment types as colData (data frame) & info on metabolites as rowData (data frame).

met.colData <- read.csv("~/Desktop/BinaryConsortiumPaper/CjEc_MetabolomicsData/MetabolomicsWorkColData.csv")
colData <- met.colData[,2:length(met.colData)]
rownames(colData) <- met.colData[,"ID"]

met.rowData <- read.csv("~/Desktop/BinaryConsortiumPaper/CjEc_MetabolomicsData/MetabolomicsWorkRowData.csv")
rowData <- met.rowData[,2:length(met.rowData)]
rownames(rowData) <- met.rowData[,"BIOCHEMICAL"]

met.data <- read.csv("~/Desktop/BinaryConsortiumPaper/CjEc_MetabolomicsData/MetabolomicsWorkAssay.csv")
assay <- met.data[,2:length(met.data)]
assay <- as.matrix(assay)
rownames(assay) <- met.rowData[,"BIOCHEMICAL"]
colnames(assay) <- met.colData[,"ID"]

#MetaboDiff function to make a Multi Assay Experiment
met <- create_mae(assay,rowData,colData)
met <- get_SMPDBanno(met, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)
      #FCreates SMPDB annotation if KEGG, HMDB and/or CHEBI IDs are available

#################################------------------#################################
#Make individual subassays for groups we're trying to compare.
 ##(1) Cj_Chi_T1 & Cj_Chi_T2

colData1 <- colData[colData$Treatment_withTime=='Cj_Chi_T2'|colData$Treatment_withTime=='Cj_Chi_T1',]
subassay1 <- assay[,colnames(assay)==colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="TP",label_colors=c("darkviolet","dodgerblue"))

met1 = knn_impute(met1,cutoff = 0.4)

met1 #19 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "TP",label_colors=c("darkviolet","dodgerblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="TP", label_colors=c("darkviolet","dodgerblue"))

#Data Analysis
met1 = diff_test(met1, group_factors = "TP")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
             group_factor="TP",
             label_colors=c("darkviolet","dodgerblue"),
             dm_cutoff=0.5,
             p_adjust = FALSE)
volcano_plot3(met1, 
             group_factor="TP",
             label_colors=c("darkviolet","dodgerblue"),
             dm_cutoff=0.5,
             p_adjust = TRUE)


################################################################################################33
##(2) Cj_NAG_T1 & Cj_NAG_T2

colData1 <- colData[colData$Treatment_withTime=='Cj_NAG_T1'|colData$Treatment_withTime=='Cj_NAG_T2',]
subassay1 <- assay[,colnames(assay)==colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="TP",label_colors=c("darkcyan","darkseagreen"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #11 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "TP",label_colors=c("darkcyan","darkseagreen"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="TP", label_colors=c("darkcyan","darkseagreen"))

#Data Analysis
met1 = diff_test(met1, group_factors = "TP")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("darkcyan","darkseagreen"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("darkcyan","darkseagreen"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(3) Ec_NAG_T1 & Ec_NAG_T2

colData1 <- colData[colData$Treatment_withTime=='Ec_NAG_T2'|colData$Treatment_withTime=='Ec_NAG_T1',]
subassay1 <- assay[,colnames(assay)==colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="TP",label_colors=c("chartreuse","goldenrod"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #15 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "TP",label_colors=c("chartreuse","goldenrod"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="TP", label_colors=c("chartreuse","goldenrod"))

#Data Analysis
met1 = diff_test(met1, group_factors = "TP")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("chartreuse","goldenrod"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("chartreuse","goldenrod"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(4) B_Chi_T1 & B_Chi_T2

colData1 <- colData[colData$Treatment_withTime=='B_Chi_T1'|colData$Treatment_withTime=='B_Chi_T2',]
subassay1 <- assay[,colnames(assay)==colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="TP",label_colors=c("yellowgreen","deepskyblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #19 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "TP",label_colors=c("yellowgreen","deepskyblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="TP", label_colors=c("yellowgreen","deepskyblue"))

#Data Analysis
met1 = diff_test(met1, group_factors = "TP")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("yellowgreen","deepskyblue"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("yellowgreen","deepskyblue"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(5) B_NAG_T1 & B_NAG_T2

colData1 <- colData[colData$Treatment_withTime=='B_NAG_T1'|colData$Treatment_withTime=='B_NAG_T2',]
subassay1 <- assay[,colnames(assay)==colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="TP",label_colors=c("yellowgreen","deepskyblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #15 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "TP",label_colors=c("yellowgreen","deepskyblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="TP", label_colors=c("yellowgreen","deepskyblue"))

#Data Analysis
met1 = diff_test(met1, group_factors = "TP")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("yellowgreen","deepskyblue"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="TP",
              label_colors=c("yellowgreen","deepskyblue"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(6) Cj_Chi_T2 & Cj_NAG_T2

colData1_1 <- colData[colData$Treatment_withTime=='Cj_Chi_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='Cj_NAG_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Chi_NAG",label_colors=c("dodgerblue","darkseagreen"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #19 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Chi_NAG",label_colors=c("dodgerblue","darkseagreen"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Chi_NAG", label_colors=c("dodgerblue","darkseagreen"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Chi_NAG")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Chi_NAG",
              label_colors=c("dodgerblue","darkseagreen"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Chi_NAG",
              label_colors=c("dodgerblue","darkseagreen"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(7) B_Chi_T2 & B_NAG_T2

colData1_1 <- colData[colData$Treatment_withTime=='B_Chi_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='B_NAG_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Chi_NAG",label_colors=c("red","deeppink"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #19 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Chi_NAG",label_colors=c("red","deeppink"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Chi_NAG", label_colors=c("red","deeppink"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Chi_NAG")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Chi_NAG",
              label_colors=c("red","deeppink"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Chi_NAG",
              label_colors=c("red","deeppink"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(8) Cj_Chi_T2 & B_Chi_T2

colData1_1 <- colData[colData$Treatment_withTime=='Cj_Chi_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='B_Chi_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("dodgerblue","red"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #18 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("dodgerblue","red"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("dodgerblue","red"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Cj_Ec")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("dodgerblue","red"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("dodgerblue","red"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(9) Cj_NAG_T2 & Ec_NAG_T2

colData1_1 <- colData[colData$Treatment_withTime=='Cj_NAG_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='Ec_NAG_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkseagreen","goldenrod"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #16 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkseagreen","goldenrod"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkseagreen","goldenrod"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Cj_Ec")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("darkseagreen","goldenrod"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("darkseagreen","goldenrod"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(10) Cj_NAG_T2 & B_NAG_T2

colData1_1 <- colData[colData$Treatment_withTime=='Cj_NAG_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='B_NAG_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkseagreen","deeppink"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #12 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkseagreen","deeppink"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkseagreen","deeppink"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Cj_Ec")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("darkseagreen","deeppink"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("darkseagreen","deeppink"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##(11) Ec_NAG_T2 & B_NAG_T2

colData1_1 <- colData[colData$Treatment_withTime=='Ec_NAG_T2',]
colData1_2 <- colData[colData$Treatment_withTime=='B_NAG_T2',]
colData1 <- rbind(colData1_1,colData1_2)
subassay1_1 <- assay[,(colnames(assay)==colData1_1$ID)]
subassay1_2 <- assay[,(colnames(assay)==colData1_2$ID)]
subassay1 <- cbind(subassay1_1,subassay1_2)

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("goldenrod","deeppink"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #15 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("goldenrod","deeppink"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("goldenrod","deeppink"))

#Data Analysis
met1 = diff_test(met1, group_factors = "Cj_Ec")
str(metadata(met1), max.level=1)

source("volcano_plot3.R")

volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("goldenrod","deeppink"),
              dm_cutoff=0.5,
              p_adjust = FALSE)
volcano_plot3(met1, 
              group_factor="Cj_Ec",
              label_colors=c("goldenrod","deeppink"),
              dm_cutoff=0.5,
              p_adjust = TRUE)


################################################################################################33
##Lookat at all T1s

colData1<- colData[colData$TP=='TP1',]
subassay1 <- assay[,colnames(assay) %in% colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #16 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkviolet","deeppink","dodgerblue"))

source("http://peterhaschke.com/Code/multiplot.R")
multiplot(
  pca_plot(met1,
           group_factor="Cj_Ec",
           label_colors=c("darkviolet","deeppink","dodgerblue")),
  tsne_plot(met1,
            group_factor="Cj_Ec",
            label_colors=c("darkviolet","deeppink","dodgerblue")),
  cols=2)

pca_plot(met1,
         group_factor="Cj_Ec",
         label_colors=c("darkviolet","deeppink","dodgerblue"))


##Add to look at NAGs too

colData1<- colData[colData$TP=='TP1',]
colData2 <- colData1[colData1$Chi_NAG=='NAG',]
subassay1 <- assay[,colnames(assay) %in% colData2$ID]

met1 <- create_mae(subassay1,rowData,colData2)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #14 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkviolet","deeppink","dodgerblue"))

source("http://peterhaschke.com/Code/multiplot.R")
multiplot(
  pca_plot(met1,
           group_factor="Cj_Ec",
           label_colors=c("darkviolet","deeppink","dodgerblue")),
  tsne_plot(met1,
            group_factor="Cj_Ec",
            label_colors=c("darkviolet","deeppink","dodgerblue")),
  cols=2)

pca_plot(met1,
         group_factor="Cj_Ec",
         label_colors=c("darkviolet","deeppink","dodgerblue"))


################################################################################################33
##Lookat at all T2s

colData1<- colData[colData$TP=='TP2',]
subassay1 <- assay[,colnames(assay) %in% colData1$ID]

met1 <- create_mae(subassay1,rowData,colData1)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #19 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkviolet","deeppink","dodgerblue"))

source("http://peterhaschke.com/Code/multiplot.R")
multiplot(
  pca_plot(met1,
           group_factor="Cj_Ec",
           label_colors=c("darkviolet","deeppink","dodgerblue")),
  tsne_plot(met1,
            group_factor="Cj_Ec",
            label_colors=c("darkviolet","deeppink","dodgerblue")),
  cols=2)

pca_plot(met1,
         group_factor="Cj_Ec",
         label_colors=c("darkviolet","deeppink","dodgerblue"))


##Add to look at NAGs too

colData1<- colData[colData$TP=='TP2',]
colData2 <- colData1[colData1$Chi_NAG=='NAG',]
subassay1 <- assay[,colnames(assay) %in% colData2$ID]

met1 <- create_mae(subassay1,rowData,colData2)
met1 <- get_SMPDBanno(met1, column_kegg_id=4,column_hmdb_id=5,column_chebi_id = NA)

#Visualizing the missing metabolite measurements across the samples
na_heatmap(met1,group_factor="Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"))

met1 = knn_impute(met1,cutoff = 0.4)
met1 #17 metabolites are excluded due to missing measurements in more than 40% of samples
met1@ExperimentList$imputed@NAMES #Finds remaining metabolites

outlier_heatmap(met1,group_factor = "Cj_Ec",label_colors=c("darkviolet","deeppink","dodgerblue"),k=2)

#Normalize
met1 <- normalize_met(met1)

quality_plot(met1, group_factor="Cj_Ec", label_colors=c("darkviolet","deeppink","dodgerblue"))

source("http://peterhaschke.com/Code/multiplot.R")
multiplot(
  pca_plot(met1,
           group_factor="Cj_Ec",
           label_colors=c("darkviolet","deeppink","dodgerblue")),
  tsne_plot(met1,
            group_factor="Cj_Ec",
            label_colors=c("darkviolet","deeppink","dodgerblue")),
  cols=2)

pca_plot(met1,
         group_factor="Cj_Ec",
         label_colors=c("darkviolet","deeppink","dodgerblue"))
