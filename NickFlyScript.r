
###Load R packages
# if you need to load the packages use install.packages("ggplot2")
##Import Step
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
#library(mctoolsr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(lubridate)
library(ggpubr)
#To find the file run file.choose() then paste it in the next line or use FlyDate=read.table(file.choose(),header=TRUE)
file.choose()
FlyData=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\data10.24.txt",header=TRUE)
FlyData


plot=ggplot(FlyData, aes(SiteLocation,Total,color=SiteLocation))+stat_summary(fun.y=mean,geom="point", size=2)+facet_wrap(~City)
plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme_set(theme_bw(base_size = 16))

#####Phyloseq import
#file.choose()
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\Table10.24.txt",header=TRUE)#For other column add Other to end of name
head(otufull)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\Meta10.24.txt",header=TRUE)
head(metadata)
taxa=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\SpeciesNames10.24.txt"))#For other column add Other to end of name
head(taxa)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
)#sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
OTU
row.names(OTU)#head(OTU) #should be 6 taxa and 40 samples
taxa_names(TAX)
TAX=tax_table(taxa)
TAX
colnames(TAX)=("Species")
sample_variables(physeq)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)
sample_data(physeq)
sample_data(physeq)$SiteLocation2=sample_data(physeq)$SiteLocation
sample_data(physeq)$SiteLocation = factor(sample_data(physeq)$SiteLocation, levels = c("North_1km","North_10km","East_1km","East_10km","South_1km","South_10km","West_1km","West_10km"))

RelativeAbu  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
Lansing=subset_samples(physeq, City == "Lansing")

plot_richness(physeq, x="DistanceFromCenter",color="DistanceFromCenter", measures=c("Shannon"))+facet_grid(~Month,scales="free_x")+geom_boxplot(aes(x=DistanceFromCenter, y=value, color=DistanceFromCenter), alpha=0.05)+
  ylab("Shannon Diversity")

plot_bar(physeq, "SiteLocation","Abundance", "Species")+xlab("Location")+ylab("Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid("Month~City",scales="free_x")#+scale_fill_manual(values=cbPalette)

ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="City",shape="City",label="SiteLocation")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = City)) #theme(legend.justification=c(1,0), legend.position=c(1,0))+


#PERMANOVA
GPdist=phyloseq::distance(physeq, "jaccard")
adonis(GPdist ~ DistanceFromCenter, as(sample_data(physeq), "data.frame"))



###RandomForest
JuneVAugust=subset_samples(physeq, Month != "July")
ForestData=JuneVAugust#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$Month)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important Taxa for classifying  samples\n by treatment")#\n in a string tells it to start a new line
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification


##########################################
###3.16.18
#########################################
file.choose()
FlyData=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\Percentages16Mar2018Other.csv",header =TRUE,check.names=FALSE)
head(FlyData)
FlyData$SiteLocationMonth=paste0(FlyData$City,FlyData$Location,FlyData$Month)
FlyData$Location = factor(FlyData$Location, levels = c("North1km","North10km","East1km","East10km","South1km","South10km","West1km","West10km"))
FlyData$Month = factor(FlyData$Month, levels = c("JUN","JUL","AUG"))
head(FlyData)
head(FlyData$Location)
plot=ggplot(FlyData, aes(Location,P.regina,color=Location,fill=Location))+stat_summary(fun.y=mean,geom="point", size=2)+facet_wrap(City~Month)+geom_bar(stat="identity")
plot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("P. regina (% of total collected)")
theme_set(theme_bw(base_size = 16))

av=aov( L.sericata~ City*Distance, data=FlyData)
summary(av)
posthoc <- TukeyHSD(av, 'City', conf.level=0.95)
posthoc

#########
#Phyloseq import
#########
otufull=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\AbundanceMatrix3.16.18All.txt",header=TRUE)#For other column add Other to end of name
#TaxTable3.16.18 for other
head(otufull)
metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\Meta3.16.18.txt",header=TRUE)
head(metadata)
taxa=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\SpeciesNames3.16All.txt"))
#taxa=as.matrix(read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Fly data\\SpeciesNames3.16Other.txt"))

head(taxa)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
#sets the plotting theme
OTU=otu_table(otufull, taxa_are_rows=TRUE)
OTU
row.names(OTU)#head(OTU) #should be 6 taxa and 40 samples
taxa_names(TAX)
TAX=tax_table(taxa)
TAX
colnames(TAX)=("Species")
sample_variables(physeq)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
taxa_names(TAX)=row.names(OTU)
physeq=phyloseq(OTU,TAX,sampdat)
sample_data(physeq)
physeq
physeq  = transform_sample_counts(physeq, function(x) x / sum(x)*100 )
physeq = filter_taxa(physeq, function(x) mean(x) > 1e-2, TRUE)
sample_data(physeq)$Location = factor(sample_data(physeq)$Location, levels = c("North1km","North10km","East1km","East10km","South1km","South10km","West1km","West10km"))
sample_data(physeq)$Month = factor(sample_data(physeq)$Month, levels = c("JUN","JUL","AUG"))

plot_richness(physeq, x="LandCover", measures=c("Shannon"))+geom_boxplot(aes(x=LandCover, y=value ), alpha=0.05)+
  ylab("Shannon Diversity")#+facet_grid(Month~City,scales="free_x") color=LandCover

plot_bar(physeq, "Location","Abundance", "Species")+xlab("Location")+ylab("Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_grid("Month~City",scales="free_x")#+scale_fill_manual(values=cbPalette)


ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="City",shape="City",label="City")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)#+facet_wrap(~City)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = City)) #theme(legend.justification=c(1,0), legend.position=c(1,0))+




df <- psmelt(physeq)
p <- ggbarplot(df, x = "City", y = "Abundance",add = c("mean_se"),#"mean_se"
               color = "black", palette = "cbPalette", facet.by="Species",
               line.color = "gray", line.size = 0.4, short.panel.labs = TRUE, p.adjust.method = "bonferroni", fill= "City") + stat_compare_means(aes(group = City), label = "..p.signif..",label.y = 7) 

p+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Relative abundance (> 1%)")+ theme(legend.position="none")

Means=compare_means(Abundance ~ City, data = df, 
                    group.by = "Species", p.adjust.method = "fdr")
#head(Means)
keeps <- c("Species","group1","group2","p.format","p.adj","method","p.signif")
keeps=Means[keeps]
#keeps


test3 <- list('Species'= keeps$Species,'group1'=keeps$group1,'group2'= keeps$group2 ,'p'=keeps$p.format,'p.adj'=keeps$p.adj,p.signif=keeps$p.signif,'Method'=keeps$method)
test3= as.data.frame(test3)
#test3
FilteredResults<-test3[!(test3$p.adj>0.2),]            
FilteredResults

#####################33
