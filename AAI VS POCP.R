##ABSOLUT PATH WHERE THE FILES ARE
path_POCP_AAI<-"C:/Users/raul/Desktop/results_POCP"


#========================================================DO NOT CHANGE ANYTHING FROM HERE=============================================
#Load packages
library(dplyr)
library(tidyr)
library(vioplot)
library(ggplot2)
library(stringr)
library(car)
library(xlsx)

#Creation of results folder
if(!dir.exists(paste0(path_POCP_AAI, "/results")))
{dir.create(paste0(path_POCP_AAI, "/results"))}


#POCP
filenames <- list.files(paste0(path_POCP_AAI, "/pocp"), pattern="*.tbl", full.names=F)

for(POCP_table in filenames){
  POCP_t<-read.delim(paste0(path_POCP_AAI,"/pocp/",POCP_table))
  n_genomes<-dim(POCP_t)[1]
  POCP_t<-gather(POCP_t, ID2, POCP_value, 2:dim(POCP_t)[2], factor_key=F)
  POCP_t<-POCP_t %>% mutate_all(list(~ str_replace(., ".faa", "")))
  colnames(POCP_t)<-c("ID1", "ID2", "POCP_value")
  POCP_t<-POCP_t[!(POCP_t$POCP_value==100 | POCP_t$POCP_value=="~"),]
  POCP_t[, 3]<- sapply(POCP_t[, 3], as.numeric)
  name<-str_remove(POCP_table, "_POCP.tbl")
  table_max_fam<-data.frame(fam_name=name, min=quantile(POCP_t$POCP_value)[1], Q2=quantile(POCP_t$POCP_value)[2], "median"=quantile(POCP_t$POCP_value)[3],Q3=quantile(POCP_t$POCP_value)[4],
                            max=quantile(POCP_t$POCP_value)[5], "mean"=mean(POCP_t$POCP_value), "SD"= sd(POCP_t$POCP_value),n_genomes=n_genomes, row.names = NULL)
  
  if(!exists("POCP_paiws"))
  {
    POCP_paiws<-POCP_t
  }else{
    POCP_paiws<-rbind(POCP_paiws,POCP_t)
  }
  if(!exists("POCP_fam"))
  {
    POCP_fam<-table_max_fam
  }else{
    POCP_fam<-rbind(POCP_fam,table_max_fam)
  }
  
}
pdf(paste0(path_POCP_AAI, "/results/POCP distribution.pdf"))
plot_dens<-density(POCP_paiws$POCP_value)
plot(plot_dens, main="POCP distribution")
polygon(plot_dens, col="darkgreen", border="blue")
vioplot(POCP_paiws$POCP_value, names=c("POCP"),
        col="gold")
max(POCP_paiws$POCP_value)
min(POCP_paiws$POCP_value)
mean(POCP_paiws$POCP_value)
sd(POCP_paiws$POCP_value)
median(POCP_paiws$POCP_value)
quantile(POCP_paiws$POCP_value)
quantile(POCP_paiws$POCP_value, probs = seq(0, 1.0, by = .01))
ggplot(data = POCP_paiws, aes(POCP_paiws$POCP_value, color = 'POCP')) +
  geom_histogram(position = "stack", bins = 100)+ labs(x = "POCP (%)")+theme_classic()
dev.off()

write.xlsx(POCP_fam,paste0(path_POCP_AAI, "/results/POCP_by_fam_GTDB_214.xlsx"), col.names = TRUE, row.names = F, append = FALSE)

#AAI
filenames_AAI <- list.files(paste0(path_POCP_AAI, "/aai"), pattern="*.tsv", full.names=F)

for(AAI_table in filenames_AAI){
  AAI_t<-read.delim(paste0(path_POCP_AAI,"/aai/",AAI_table))
  AAI_t<-subset(AAI_t, select = c("Label.1", "Label.2", "AAI", "Proteome.cov."))
  colnames(AAI_t)<-c("ID1", "ID2", "AAI_value", "Coverage_AAI")
  matrix_AAI<- pivot_wider(subset(AAI_t, select = c("ID1", "ID2", "AAI_value")), names_from = "ID2", values_from = "AAI_value")
  matrix_coverage<-pivot_wider(subset(AAI_t, select = c("ID1", "ID2", "Coverage_AAI")), names_from = "ID2", values_from = "Coverage_AAI")
  first_col_order <- order(matrix_AAI$ID1)
  matrix_AAI <- data.frame(matrix_AAI[, c(1, first_col_order + 1)])
  matrix_coverage<-data.frame(matrix_coverage[, c(1, first_col_order + 1)])
  matrix_AAI[,2:ncol(matrix_AAI)][upper.tri(matrix_AAI[,2:ncol(matrix_AAI)])] <- NA
  matrix_coverage[,2:ncol(matrix_coverage)][upper.tri(matrix_coverage[,2:ncol(matrix_coverage)])] <- NA
  AAI_t <- matrix_AAI %>% gather("ID2", "AAI_value", -ID1)
  AAI_t$AAI_value<-as.numeric(AAI_t$AAI_value)
  AAI_t$Coverage_AAI<-as.numeric((matrix_coverage%>% gather("ID2", "Coverage_AAI", -ID1))[,3])
  AAI_t<-AAI_t[!(AAI_t$AAI_value==100.0|is.na(AAI_t$AAI_value)),]
  n_genomes<-length(unique(c(AAI_t$ID1,AAI_t$ID2)))
  name<-str_remove_all(AAI_table, ".tsv|aai_")
  table_max_fam<-data.frame(fam_name=name, min=quantile(AAI_t$AAI_value)[1], Q2=quantile(AAI_t$AAI_value)[2], "median"=quantile(AAI_t$AAI_value)[3],Q3=quantile(AAI_t$AAI_value)[4],
                            max=quantile(AAI_t$AAI_value)[5], "mean"=mean(AAI_t$AAI_value), "SD"= sd(AAI_t$AAI_value),n_genomes=n_genomes, row.names = NULL)
  
  if(!exists("AAI_paiws"))
  {
    AAI_paiws<-AAI_t
  }else{
    AAI_paiws<-rbind(AAI_paiws,AAI_t)
  }
  if(!exists("AAI_fam"))
  {
    AAI_fam<-table_max_fam
  }else{
    AAI_fam<-rbind(AAI_fam,table_max_fam)
  }
  
}

plot_dens<-density(AAI_paiws$AAI_value)
plot_cov <- density(AAI_paiws$Coverage_AAI)
pdf(paste0(path_POCP_AAI, "/results/AAI distribution.pdf"))
plot(plot_dens, main="AAI distribution")
polygon(plot_dens, col="darkgreen", border="blue")
plot(plot_cov, main="Proteome Coverage")
polygon(plot_cov, col="dark red", border="blue")
vioplot(AAI_paiws$AAI_value, names=c("AAI"),
        col="orange")
max(AAI_paiws$AAI_value)
min(AAI_paiws$AAI_value)
mean(AAI_paiws$AAI_value)
sd(AAI_paiws$AAI_value)
median(AAI_paiws$AAI_value)
quantile(AAI_paiws$AAI_value)
quantile(AAI_paiws$AAI_value, probs = seq(0, 1.0, by = .01))
ggplot(data = AAI_paiws, aes(AAI_paiws$AAI_value, color = 'AAI')) +
  geom_histogram(position = "stack", bins = 100)+ labs(x = "AAI (%)")+theme_classic()
dev.off()
#Relatioship between AAI and the COVERAGE
pdf(paste0(path_POCP_AAI, "/results/AAI VS coverage.pdf"))
scatterplot(AAI_value ~ Coverage_AAI, data = AAI_paiws, col= "darkolivegreen4",
            xlab="Coverage (%)", ylab="AAI (%)", regLine = list(col="darkgreen"), smooth=list(col.smooth="orange", col.spread="orange"))
dev.off()
regresion <- lm(AAI_value ~ Coverage_AAI, data = AAI_paiws)
summary(regresion)
max(AAI_paiws$Coverage_AAI)
min(AAI_paiws$Coverage_AAI)
mean(AAI_paiws$Coverage_AAI)
median(AAI_paiws$Coverage_AAI)
sd(AAI_paiws$Coverage_AAI)

write.xlsx(AAI_fam,paste0(path_POCP_AAI, "/results/AAI_by_fam_GTDB_214.xlsx"), col.names = TRUE, row.names = F, append = FALSE)


#comparison of the two
AAI_POCP_table<-merge(POCP_paiws, AAI_paiws, by.x = c("ID1", "ID2"), by.y = c("ID2", "ID1"), all.x = T)
dens_comp <- gather(AAI_POCP_table, key = 'analysis', value = 'value', contains('value'))
dens_comp$analysis <-str_remove(dens_comp$analysis, "_value")
pdf(paste0(path_POCP_AAI, "/results/POCP VS distribution.pdf"))
ggplot(dens_comp, aes(x=value, fill=analysis)) +
  geom_density(aes(y = after_stat(count)), alpha=0.25)+ scale_fill_manual(values=c("red4", "darkgreen"))+ theme_classic() 
dev.off()

regresion <- lm(POCP_value ~ AAI_value, data = AAI_POCP_table)
summary(regresion)
pdf(paste0(path_POCP_AAI, "/results/AAI VS POCP.pdf"))
scatterplot(AAI_value ~ POCP_value, data = AAI_POCP_table, col= "royalblue3",
            xlab="POCP (%)", ylab="AAI (%)", regLine = list(col="red4"), smooth=list(col.smooth="orange", col.spread="orange"))
dev.off()
write.xlsx(AAI_POCP_table,paste0(path_POCP_AAI, "/results/POCP_vs_AAI_pairwise_GTDB_214.xlsx"), col.names = TRUE, row.names = F, append = FALSE)

