install.packages("vcfR")
install.packages("adegenet")
library(vcfR)
library(adegenet)


###importing file to vcfR and creating genid 
vcf <- read.vcfR("tilapia.vcf", verbose=FALSE)
show(vcf)
gid <- vcfR2genind(vcf,return.alleles=TRUE)
gid
save(gid, file="gid.RData")

#/// GENIND OBJECT /////////
  
#  // 379 individuals; 60,785 loci; 111,871 alleles; size: 192.7 Mb

#// Basic content
#@tab:  379 x 111871 matrix of allele counts
#@loc.n.all: number of alleles per locus (range: 1-2)
#@loc.fac: locus factor for the 111871 columns of @tab
#@all.names: list of allele names for each locus
#@ploidy: ploidy of each individual  (range: 2-2)
#@type:  codom
#@call: adegenet::df2genind(X = t(x), sep = sep)

#// Optional content
#- empty -

  
## Load packages
install.packages("pegas")
install.packages("poppr")
install.packages("miscTools")
install.packages("FactoMineR")
lib = c("adegenet","ggplot2","RColorBrewer","pegas","poppr","miscTools","stringr","FactoMineR")
lapply(lib, library, character.only=TRUE)

###############################################
##############################################
###remove non-polymorphic markers

load("gid.RData")
gid
summary(gid$pop)
indNames(gid)

#no names
#name your populations
pop_names <- indNames(gid)
pop2 <- rep (NA, length(pop_names))

pop2[grep("Kange", pop_names)] <- "kange"
pop2[grep("Kikota", pop_names)] <- "kikota"
pop2[grep("Victoria", pop_names)] <- "victoria"
pop2[grep("Mukoke", pop_names)] <- "mukoke"
pop2[grep("Pearl", pop_names)] <- "pearl"
pop2[grep("AA", pop_names)] <- "aa_fisheries"
pop2[grep("Rubona", pop_names)] <- "rubona"
pop2[grep("Albert", pop_names)] <- "albert"
pop2[grep("Kumi", pop_names)] <- "kumi"
pop2[grep("Blue_Valley", pop_names)] <- "blue_valley"
pop2[grep("Kyoga", pop_names)] <- "kyoga"
pop2[grep("Nile", pop_names)] <- "nile"
pop2[grep("Kasolwe", pop_names)] <- "kasolwe"
pop2[grep("Asi", pop_names)] <- "asi"
pop2


pop(gid) 
#nothing in the pop field at the moment for this genind...
# try attaching this vector of popdata to the genind object...
#?strata
strata(gid) <- data.frame(pop2)
strata(gid)
setPop(gid) <- ~pop2
gid

# so there is now pop data attached...
pop(gid)

# rename the dataset, and then the population sample names 
data_filt <- gid

data_filt
nPop(data_filt)
# shows 14 thats correct

summary(data_filt$pop)
# Pop names are fine, but use below if wanting to rename anything

data_filt
#/// GENIND OBJECT /////////
#  
#  // 379 individuals; 60,785 loci; 111,871 alleles; size: 192.7 Mb
#
#// Basic content
#@tab:  379 x 111871 matrix of allele counts
#@loc.n.all: number of alleles per locus (range: 1-2)
#@loc.fac: locus factor for the 111871 columns of @tab
#@all.names: list of allele names for each locus
#@ploidy: ploidy of each individual  (range: 2-2)
#@type:  codom
#@call: adegenet::df2genind(X = t(x), sep = sep)

#// Optional content
#@pop: population of each individual (group size range: 1-35)
#@strata: a data frame with 1 columns ( pop2 )

## Check loci for missing data
locmiss = propTyped(data_filt, by="loc")
head(locmiss)

# explore which loci meet various missing data thresholds...
locmiss[which(locmiss < 0.95)] # more than 5% missing data 
# many missing >10K

locmiss[which(locmiss < 0.951)] # more than 4.9% missing data 
# several hundred hits, so that shows that the dataset is all - 0.95 

#
data_filt = missingno(data_filt, type="loci", cutoff=0.05)

#save environment as locmiss takes ages to run
save.image(file='tilapia.RData')

# Now assess missing data by individual...

indmiss = propTyped(data_filt, by="ind")
head(indmiss)

indmiss[which(indmiss < 0.95)] 
#  1 individual missing
#Kikota_Integrated_fish_farm_18_246373_C23.CEL 
#0.9497917 

indmiss[which(indmiss < 0.9)] 
# 0

## Poppr VERY USEFUL FUNCTIONS ###

# create a new dataset, filtering out individuals which fail the missing data theshold
data_filt = missingno(data_filt, type="geno", cutoff=0.05)
## 1 individual removed

data_filt
summary(data_filt$pop)
nPop(data_filt)
# full filtered dataset is now 397 samples typed at at least 90% of 20,938 loci, across 37 pops

#check for polymorphic loci (some may be monomorphic after filtering individuals)

## Are loci polymorphic?
#?isPoly
isPoly(data_filt)
poly_snps <- isPoly(data_filt) 
summary(poly_snps)

#Mode   FALSE    TRUE 
#logical   21425   29224 

# 21425 K snps are monomorphic
# save a genind of these 397 individuals across all 20k loci (including potential monomorphs)

## Output filtered genind object 
tilapiafulldataname = paste("tilapia_",nInd(data_filt),"I-",nLoc(data_filt),
                           "L-",nPop(data_filt),"P.RData",sep="")
tilapiafulldataname
save(data_filt, file=tilapiafulldataname)
# saves genind object and shows number of samples / snps / pops

# Now create another genind with the potentially monomorphic loci removed...

# create a list of these loci which ARE polymorphic...
poly_loci = names(which(isPoly(data_filt) == TRUE))
mono_loci = names(which(isPoly(data_filt) == FALSE))
poly_loci
#mono_loci
data_filt
#250,649 loci before filtering out monomorphic loci...

# filter the dataset to retain only the polymorphic loci
data_filt = data_filt[loc=poly_loci]
data_filt
#29,224 loci, all polymorphic

# save this genind dataset

## Output filtered genind object 
tilapiafulldataname2 = paste("tilapia_poly_",nInd(data_filt),"I-",nLoc(data_filt),
                            "L-",nPop(data_filt),"P.RData",sep="")
tilapiafulldataname2
save(data_filt, file=tilapiafulldataname2)

mono_loci <- as.data.frame(mono_loci)
mono_loci
mono_name = paste(mono_loci, sep="\,")
mono_name


x <- as.data.frame(mono_loci)
x2 <- write.table(x,  sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.csv(x2, file="mono_loci.txt")
# saves genind object and shows number of samples / snps / pops




################################################################
###############DAPC analysis
################################################################
load("tilapia_poly_378I-29224L-14P.RData")
# this is the polymorphic dataset, with only 16k loci...

data_filt
# 397 inds, 37pops, 16.661k loci 

summary(data_filt$pop)
indNames(data_filt)

## Cross validation to find the optimal number of PCs to retain in dapc
# trial run with only 5 reps - re-run with 100+ reps for paper!!!

####
##as suggested by charlie
## Scale of above xval validation has optimal PCs in massuve intervals of 50 
## We can re-run with a narrower n.pca.max window & more reps to assess it at a finer scale
## AGAIN - trial below uses only 5 reps so increase to 100+ for full analysis runs...

#will use 13 (k-1 n of pops =14)

xtil = tab(data_filt, NA.method="mean")
xvaltil = xvalDapc(xtil, data_filt$pop, n.pca.max=13, training.set=0.9,
                   result="groupMean", center=TRUE, scale=FALSE,
                   n.rep=100, n.pca=NULL, parallel="snow")

## Number of PCs with best stats
xvaltil$`Number of PCs Achieving Highest Mean Success`
xvaltil$`Number of PCs Achieving Lowest MSE`
xvaltil$`Root Mean Squared Error by Number of PCs of PCA` # lower score = better
# 11 PCs to retain


## Number of PCs with best stats
xval$`Number of PCs Achieving Highest Mean Success`
#11
xval$`Number of PCs Achieving Lowest MSE`
#11
xval$`Root Mean Squared Error by Number of PCs of PCA` # lower score = better
#11


## Run the DAPC using population IDs as priors, using 80 PCs to retain
mydapc = dapc(data_filt, data_filt$pop, n.pca=11, n.da=3)


## Analyse how much percent of genetic variance is explained by each axis
percent = mydapc$eig/sum(mydapc$eig)*100
percent
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", 
        names.arg=round(percent, 2))
# good distribution in power pc1 pc2 drpo to pc3 pc4
# 38.65% of variation explained by PC1, ~28% by PC2 but only ~10% by PC3.  
# P2 & PC3 can both be plot against PC1 - but not much power expected beyond PC1v2

## SNP contributions
contrib = loadingplot(mydapc$var.contr)
contrib = loadingplot(mydapc$var.contr, threshold=0.0005)
contrib$var.names
# 0 >0.0005 loading threshold - these are the most powerful rangewide snps via DAPC 

## Create ATLdataframe with PC info
df_pca = as.data.frame(mydapc$ind.coord)

## Add a column containing individuals
df_pca = cbind(rownames(df_pca), df_pca)

## Rename first three columns
colnames(df_pca) = c("Indiv","PC1","PC2","PC3")

## Add a column with the population IDs
df_pca$pop = data_filt$pop
summary(df_pca$pop)
df_pca$Indiv
## shows pop names attributed to df_pca

## Flip axis by multiplying by minus 1 (use -2 or -3 for further rotation!)
#df_pca[ , 2:4] = df_pca[ , 2:4]*-4

# ----------------- #
#
# Visualise DAPC using all samples
#
# ----------------- #

## set ggplot2 themes

load("tilapia.RData")

ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=15),
                axis.text.x = element_text(colour="black", size=15),
                axis.title = element_text(colour="black", size=15),
                # legend.justification = c(1,1),
                legend.position = "top",
                legend.text = element_text(size=10),
                legend.key.size = unit(0.7,"cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.grid = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(linewidth = 0.2, colour="black", fill=NA),
                plot.title = element_text(hjust=0.5, size=20) # title centered 
)

# similar theme with alternative sizing for multiplot figure (for smaller publication plots)

## ggplot2 theme2
#ggtheme2 = theme(legend.title = element_blank(),
#                 axis.text.y = element_text(colour="black", size=10),
#                 axis.text.x = element_text(colour="black", size=10),
#                 axis.title = element_text(colour="black", size=15),
#                 # legend.justification = c(1,1),
#                 legend.position = "top",
#                 legend.text = element_text(size=12),
#                 legend.key.size = unit(0.7,"cm"),
#                 legend.box.spacing = unit(0, "cm"),
#                 panel.grid = element_blank(),
#                 panel.background = element_blank(),
#                 panel.border = element_rect(colour="black", fill=NA, size=1),
#                 plot.title = element_text(hjust=0.5, size=25) # title centered 
#)

## Function to add region or species labels to population labels
summary(data_filt$pop)

##make it by lake
# If pop label is present then this function will output the region/species on the plot key

addregionfull.K11C5 = function(x){
  if(x=="victoria") y = "Lake Victoria"
  if(x=="albert") y = "Lake Albert" 
  if(x=="kyoga") y = "Lake Kyoga"
  if(x=="nile") y = "Nile river"
  if(x=="kange"|x=="kikota" |x=="asi"|x=="aa_fisheries" |x=="rubona"|x=="kumi" |x=="blue_valley" |x=="kasolwe"|x=="mukoke" |x=="pearl")
    y = "Fish farms"
  return(y)
}

addregionfull.K11C5        

#addregionfull.K5 = function(x){
#  # If pop label is present function will output the country
#  if(x=="Mac"|x=="Aeg") y = " Aegean "
#  if(x=="Adr") y = " Adriatic "
#  if(x=="Csa"|x=="Laz"|x=="Sar") y = " Western Med. "
#  if(x=="Tan"|x=="Pen"|x=="Vig"|x=="Idr"|x=="Bre"|
#     x=="Jer"|x=="Loo"|x=="Sbs"|x=="Ios"|x=="Don"|
#     x=="Kil"|x=="Cor"|x=="Iom"|x=="Pem"|x=="Pad"|
#     x=="Heb"|x=="Ork"|x=="She"|x=="Brd"|x=="Cro"|
#     x=="Eye"|x=="Mul"|x=="Ven") y = " Atlantic "
#  if(x=="Hel"|x=="Flo"|x=="Lys"|x=="Ber"|x=="Tro"|x=="Sin") y = " Scandinavia "
#  return(y)
#}

df_pca$region = sapply(df_pca$pop, addregionfull.K11C5) 
unique(df_pca$region)

## Reorder levels for plotting
df_pca$region = factor(df_pca$region, 
                       levels=c( "Lake Victoria",
                                 "Lake Albert", 
                                 "Lake Kyoga",
                                 "Nile river",
                                 "Fish farms"
                       ))

## Calculate centroid position for each population
centroid = aggregate(cbind(PC1, PC2, PC3) ~ pop, data=df_pca, FUN=mean)
centroid$region = sapply(centroid$pop, addregionfull.K11C5)

## Find and store coordinate info required to draw segments
segs = merge(df_pca, setNames(centroid, 
                              c("pop","oPC1S1","oPC2S2","oPC3S3")),
             by = "pop", sort = FALSE)

## Definitions
n_snp = nLoc(data_filt)
n_ind = nInd(data_filt)
n_pop = nPop(data_filt)

## Colour Brewer definitions (n colours must match n regions/pops)
colsfullK11.C5 = c("yellow","green3","purple","red", "grey")
#colsfullK5.C14 = c("yellow","olivedrab1","green3","blue","purple","red","orange","hotpink","cyan",
#                  "coral","chocolate","chartreuse","brown4","bisque")


# PLOT DAPCS!
# First up, PC1 vs PC2

## Scatter plot axis 1 vs 2
DAPCfull_PC12=ggplot(df_pca, aes(x=PC1, y=PC2))+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spiders
  geom_segment(data=segs, mapping = aes(xend=oPC1S1, yend=oPC2S2,
                                        colour=region), show.legend=FALSE)+
  geom_point(aes(fill=region), shape=21, size=2, show.legend=TRUE)+
  # centroids
  geom_label(data=centroid, size=3, aes(label=pop, fill=region), show.legend=FALSE, position=position_jitter())+
  #geom_text(size = 3, aes (label=pop))
  scale_fill_manual(values=colsfullK11.C5)+
  scale_color_manual(values=colsfullK11.C5)+
  labs(x=paste("PC1 (",round(percent[1],digits=1)," %)", sep=""))+
  labs(y=paste("PC2 (",round(percent[2],digits=1)," %)", sep=""))+
  # ggtitle(paste("European lobster DAPC:", no_ind,"individuals |",
  #               no_snps, "SNPs"))+
  ggtheme+
  # labs(tag="A")+
  # legend key on two lines
  guides(fill=guide_legend(nrow=1, byrow=TRUE))

DAPCfull_PC12

## 

## Export plot
titlefulldapc = paste("dapc_fullK11.C5_",n_ind,"I-",
                      n_snp,"L-",n_pop,"P", sep="")
titlefulldapc
ggsave(paste(titlefulldapc,"_PC1-2.png",sep=""), width=12, height=12, dpi=300)
#ggsave(paste(titleATLfulldapc,"PC1-2_K3_.pdf",sep=""), width=12, height=12)

scatter(mydapc)
## Now Scatter plot axis 1 vs 3
DAPCfull_PC13=ggplot(df_pca, aes(x=PC1, y=PC3))+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spiders
  geom_segment(data=segs, mapping = aes(xend=oPC1S1, yend=oPC3S3,
                                        colour=pop), show.legend=FALSE)+
  geom_point(aes(fill=pop), shape=21, size=3, show.legend=TRUE)+
  # centroids
  geom_label(data=centroid, size=5, aes(label=pop, fill=pop), show.legend=FALSE)+
  scale_fill_manual(values=colsfullK5.C9)+
  scale_color_manual(values=colsfullK5.C9)+
  labs(x=paste("PC1 (",round(percent[1],digits=1)," %)", sep=""))+
  labs(y=paste("PC3 (",round(percent[3],digits=1)," %)", sep=""))+
  # ggtitle(paste("European lobster DAPC:", no_ind,"individuals |",
  #               no_snps, "SNPs"))+
  ggtheme+
  # labs(tag="A")+
  # legend key on two lines
  guides(fill=guide_legend(nrow=3, byrow=TRUE))

DAPCfull_PC13

## 

## Export plot
ggsave(paste(titlefulldapc,"_PC1-3.png",sep=""), width=12, height=12, dpi=300)
#ggsave(paste(titleneutraldapc,"PC1-3_K5_.pdf",sep=""), width=12, height=12)


#low tide vs high tide

# create secondary dataset minus the other outliers...
summary(data_filt$pop)
data_sub = popsub(data_filt, sublist=c("H"))
data_sub
nPop(data_sub)
summary(data_sub$pop)
#rerun the whole analysis


  