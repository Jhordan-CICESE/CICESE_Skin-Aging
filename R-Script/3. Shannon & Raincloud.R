library(matrixStats)
library(pipeR)
library(rlist)
library(purrr)
library(ggpubr)
library(dplyr)
library(plotly)
library(vegan) # Shannon
library(ggplot2) ## plotting
library(systemfonts) ## custom fonts
library(lawstat) # p-value Shannon
library(outliers) # p-value Shannon

## INSTALL PACKAGES ----------------------------------------------
## install CRAN packages if needed
pckgs <- c("ggplot2", "systemfonts", "ggforce", 
           "ggdist", "ggbeeswarm", "devtools")
new_pckgs <- pckgs[!(pckgs %in% installed.packages()[,"Package"])]
if(length(new_pckgs)) install.packages(new_pckgs)

## install gghalves from GitHub if needed
if(!require(gghalves)) {  
  devtools::install_github('erocoar/gghalves')
}

## CUSTOM THEME --------------------------------------------------
## overwrite default ggplot2 theme
theme_set(
  theme_minimal(
    ## increase size of all text elements
    base_size = 18, 
    ## set custom font family for all text elements
    base_family = "Oswald")
)

## overwrite other defaults of theme_minimal()
theme_update(
  ## remove major horizontal grid lines
  panel.grid.major.x = element_blank(),
  ## remove all minor grid lines
  panel.grid.minor = element_blank(),
  ## remove axis titles
  axis.title = element_blank(),
  ## larger axis text for x
  axis.text.x = element_text(size = 16),
  ## add some white space around the plot
  plot.margin = margin(rep(8, 4))
)

WD <- "G:\\Mi unidad\\CICESE\\Tesis\\2. Analysis\\2. Results\\2. FBMN\\2. ChemicalClasses\\2. RainCloud"
setwd(WD)

# This table was manually revised to add Pubchem ID
Final_Table <- read.csv(paste0(WD, "//0. Source//1. Final-Annotation-Enriched_FBMN.csv"))

#read MS1 data
ms1_data_normalyzed <- read.csv(paste0(WD, "//0. Source//2. Data_Norm_Imp_Log.csv"))

#To perform raincloud plots
met_complete <- merge(ms1_data_normalyzed, Final_Table,
    by.x="Feature", by.y="Feature")

met_complete_sub <- subset(met_complete, 
    select = c(1:76,83:85))

met_charts <- data.frame(t(met_complete_sub[,-1]))
#colnames(met_charts) <- met_complete_sub$Feature

met_charts$Cohort  <- c("Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Jóvenes","Mayores","Jóvenes",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores","Mayores","Mayores","Mayores",
                        "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                        "Jóvenes","Jóvenes","Mayores","Jóvenes","Jóvenes",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores","Mayores","Mayores","Mayores","Mayores",
                        "Mayores",
                        "Subclass", "Class", "Superclass")

#Filter by cohort
met_charts_Young <- met_charts %>% dplyr::filter( Cohort == "Jóvenes")
met_charts_Old <- met_charts %>% dplyr::filter( Cohort == "Mayores")

#Extract metabolites class
met_class <- met_charts[77,]

#Bind it to the subgrouped
met_charts_Young <- bind_rows(met_charts_Young, met_class)
met_charts_Old <- bind_rows(met_charts_Old, met_class)

#Transpose to make pie charts
met_charts_Young <- data.frame(t(met_charts_Young))
met_charts_Old <- data.frame(t(met_charts_Old))

#Remove row names
rownames(met_charts_Young) <- NULL
rownames(met_charts_Old) <- NULL

#Convert to numeric
met_charts_Young <- mutate_all(met_charts_Young, function(x) as.numeric(as.character(x)))
met_charts_Old <- mutate_all(met_charts_Old, function(x) as.numeric(as.character(x)))

#transpose metabolite classe to re insert it
met_class <-t(met_class)

met_charts_Young$Class <- c(met_class)
met_charts_Old$Class <- c(met_class)

#Remove last rows to don't have NA after numeric convertion
met_charts_Young <- met_charts_Young[-c(2600),]
met_charts_Old <- met_charts_Old[-c(2600),]

#Row Sum
met_charts_Young$Total <- rowSums(met_charts_Young[,1:38])
met_charts_Old$Total <- rowSums(met_charts_Old[,1:37])


#Group by metabolite class
met_Young_statistics <- aggregate(met_charts_Young[,1:38], by=list(Class=met_charts_Young$Class), FUN=sum)
met_Old_statistics <- aggregate(met_charts_Old[,1:37], by=list(Class=met_charts_Old$Class), FUN=sum)

#Row Mean
met_Young_statistics$Mean <- rowMeans(met_Young_statistics[,2:39])
met_Old_statistics$Mean <- rowMeans(met_Old_statistics[,2:38])

#Row SD
met_Young_statistics$SD <- rowSds(as.matrix(met_Young_statistics[,2:39]))
met_Old_statistics$SD <- rowSds(as.matrix(met_Old_statistics[,2:38]))

met_statistics_Cohort <- full_join(met_Young_statistics, met_Old_statistics, by = "Class")

T_met_statistics_Cohort <- data.frame(t(met_statistics_Cohort))
colnames(T_met_statistics_Cohort) <- met_statistics_Cohort$Class
T_met_statistics_Cohort <- T_met_statistics_Cohort[-c(1,40,41,79,80),]

#Convert to numeric
T_met_statistics_Cohort <- mutate_all(T_met_statistics_Cohort, function(x) as.numeric(as.character(x)))
T_met_statistics_Cohort[is.na(T_met_statistics_Cohort)] <- 0

Cohorte <- c("Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Jóvenes","Jóvenes","Jóvenes","Jóvenes",
            "Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores","Mayores","Mayores","Mayores",
            "Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores","Mayores","Mayores","Mayores","Mayores",
            "Mayores")

#T_met_statistics_Cohort <- T_met_statistics_Cohort[,-c(1,2)]

T_met_statistics_Cohort <-cbind(Cohorte, T_met_statistics_Cohort)

test.fun <- function(T_met_statistics_Cohort, col) { 

 c1 <- combn(unique(T_met_statistics_Cohort$Cohort),2)
 sigs <- list()
 for(i in 1:ncol(c1)) {
    sigs[[i]] <- t.test(
                   T_met_statistics_Cohort[T_met_statistics_Cohort$Cohort == c1[1,i],col],
                   T_met_statistics_Cohort[T_met_statistics_Cohort$Cohort == c1[2,i],col]
                 )
    }
    names(sigs) <- paste("Cohorte",c1[1,],"vs Cohorte",c1[2,])

 tests_Cohort <- data.frame(Test=names(sigs),
                    W=unlist(lapply(sigs,function(x) x$statistic)),
                    p=unlist(lapply(sigs,function(x) x$p.value)),row.names=NULL)

 return(tests_Cohort)
}

tests_Cohort <- lapply(colnames(T_met_statistics_Cohort)[-c(1,1)],function(x) test.fun(T_met_statistics_Cohort,x))

names(tests_Cohort) <- paste0((colnames(T_met_statistics_Cohort)[-c(1,1)]),"_")

significative_class_Young_vs_Old <- map(tests_Cohort, ~ dplyr::filter(.x,
Test == "Jóvenes vs Mayores"))

significative_class_Young_vs_Old <- significative_class_Young_vs_Old %>>% 
list.filter(p < 0.05)

write.csv(significative_class_Young_vs_Old,".//3. Statistics//result_T-test_Young_vs_Old.csv",
          row.names = FALSE)
write.csv(T_met_statistics_Cohort,".//3. Statistics//Abundance per class_Cohort.csv",
          row.names = FALSE)

#Convert to numeric
met_comp_sub <- mutate_all(met_complete_sub, function(x) as.numeric(as.character(x)))

met_comp_sub$X <- met_complete_sub$X
met_comp_sub$Class <- met_complete_sub$Class
met_comp_sub$Subclass <- met_complete_sub$Subclass
met_comp_sub$Superclass <- met_complete_sub$Superclass
met_comp_sub$Feature <- NULL
met_comp_sub <- met_comp_sub[,-c(76,78)]

#Row Sum
met_comp_sub$Total_Young <- rowSums( met_comp_sub[,c(1:18,20,40:56,58,59)])
met_comp_sub$Total_Old <- rowSums( met_comp_sub[,c(19,21:39,57,60:75)])

#Row Mean
met_comp_sub$Mean_Young <- rowMeans( met_comp_sub[,c(1:18,20,40:56,58,59)])
met_comp_sub$Mean_Old <- rowMeans( met_comp_sub[,c(19,21:39,57,60:75)])

#Row SD
met_comp_sub$SD_Young <- rowSds(as.matrix( met_comp_sub[,c(1:18,20,40:56,58,59)]))
met_comp_sub$SD_Old <- rowSds(as.matrix( met_comp_sub[,c(19,21:39,57,60:75)]))

#Classes with FC>1.5 piechart
FC2_classes <- data.frame(class = c(met_comp_sub[,76]),
                          Mean_Young = c(met_comp_sub[,79]), Mean_Old = c(met_comp_sub[,80]))

FC2_classes$FC_YoungvOld <- (FC2_classes$Mean_Young/FC2_classes$Mean_Old)
FC2_Young_vs_Old <- FC2_classes %>% dplyr::filter(FC_YoungvOld >= 1.5)
FC2_Young_vs_Old <- FC2_Young_vs_Old[,-c(4,6)]

#Plotly Module----
id <- "1sZFtGZfnqGzvAhUN32j0pg614MaEm3lN"
class_colors <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", id), na.strings=c("", "NA"))

Classes_FC_Young_vs_Old <- FC2_Young_vs_Old %>% dplyr::count(class)
Classes_FC_Young_vs_Old <- merge(Classes_FC_Young_vs_Old,class_colors, by.x = "class", by.y = "X.canopus.Class.") 


fig_class_Young_vs_Old <- plot_ly(Classes_FC_Young_vs_Old, labels = ~class, values = ~n, type = 'pie', 
                                  sort = FALSE, 
                                  direction = "clockwise",
                                  marker = list(colors= ~mycolors441), textfont = list(size = 15))%>%
  layout(font=list(size = 12),
         legend = list(title=list(text='<b> Class </b>'), x= 100, y=0.5))
fig_class_Young_vs_Old

#Set Enviroments
Sys.setenv("plotly_username" = "froz9")
Sys.setenv("plotly_api_key" = "3YcsdCmBtC8rfYLZA3PR")

#Save plots
plotly_IMAGE(fig_class_Young_vs_Old, format="png", width = 1000, height = 1000, out_file=".//2. Plots//Class_FC_1.5_Young_vs_Old.png")

write.csv(FC2_Young_vs_Old,paste0(".//3. Statistics//FC_1.5_Young_vs_Old.csv"),
          row.names = FALSE)
write.csv(Classes_FC_Young_vs_Old,paste0(".//3. Statistics//Classes_Young_vs_Old.csv"),
          row.names = FALSE)

#--------------- CALCULATE SHANNON------------------------------#
sha <- T_met_statistics_Cohort[,-1]

#Calculate the diversity indexes
Diversity_Shannon <- diversity(sha,index = "shannon")

Results <- data.frame("Shannon" = unlist(Diversity_Shannon))

Results$Muestra <- rownames(T_met_statistics_Cohort)

#By Sexo
Results$Sexo  <- c("Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres",
                     "Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Hombres",
                     "Hombres","Hombres","Hombres","Hombres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres","Mujeres","Mujeres","Mujeres","Mujeres",
                     "Mujeres")

Boxplot_Shannon <- ggplot(Results, aes(x = Sexo, y = Shannon,
                                       fill= Sexo)) + geom_violin() + geom_point(position = position_jitter(seed = 1, width = 0.2))
Boxplot_Shannon

ggsave(".//1. Shannon//Shannon_boxplot_Sexo.png",
       device = "png",
       width = 2500,
       height = 1250,
       units = "px",
       dpi = "retina")

#By Fourties
Results$Edad  <- c("Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Jóvenes","Jóvenes",
                       "Jóvenes","Jóvenes","Jóvenes","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores","Mayores","Mayores","Mayores","Mayores",
                       "Mayores")
Boxplot_Shannon <- ggplot(Results, aes(x = Edad, y = Shannon,
                                       fill= Edad)) + geom_violin() + geom_point(position = position_jitter(seed = 1, width = 0.2))
Boxplot_Shannon

ggsave(".//1. Shannon//2. Fourties//Shannon_boxplot_Edad.png",
       device = "png",
       width = 2500,
       height = 1250,
       units = "px",
       dpi = "retina")

#By Cohorte
Results$Cohorte  <- c("H-J","H-J","H-J","H-J","H-J",
                    "H-J","H-J","H-J","H-J","H-J",
                    "H-J","H-J","H-J","H-J","H-J",
                    "H-J","H-J","H-J","H-J",
                    "M-J","M-J","M-J","M-J","M-J",
                    "M-J","M-J","M-J","M-J","M-J",
                    "M-J","M-J","M-J","M-J","M-J",
                    "M-J","M-J","M-J","M-J",
                    "H-M","H-M","H-M","H-M","H-M","H-M",
                    "H-M","H-M","H-M","H-M","H-M",
                    "H-M","H-M","H-M","H-M","H-M",
                    "H-M","H-M","H-M","H-M",
                    "M-M","M-M","M-M","M-M","M-M","M-M",
                    "M-M","M-M","M-M","M-M","M-M",
                    "M-M","M-M","M-M","M-M","M-M",
                    "M-M")

Boxplot_Shannon <- ggplot(Results, aes(x = Cohorte, y = Shannon,
                                       fill= Cohorte)) + geom_violin() + geom_point(position = position_jitter(seed = 1, width = 0.2))
Boxplot_Shannon

ggsave(".//1. Shannon//2. Fourties//Shannon_boxplot_Cohorte.png",
       device = "png",
       width = 2500,
       height = 1250,
       units = "px",
       dpi = "retina")



#ANALISIS ESTADISTICO
summary(Results$Shannon)
#Variables
Hombres <- Results %>% filter(Results$Sexo == "Hombres") %>% pull(Shannon)
Mujeres <- Results %>% filter(Results$Sexo == "Mujeres") %>% pull(Shannon)

Jovenes <- Results %>% filter(Results$Edad == "Jóvenes") %>% pull(Shannon)
Mayores <- Results %>% filter(Results$Edad == "Mayores") %>% pull(Shannon)

HJ <- Results %>% filter(Results$Cohorte == "H-J") %>% pull(Shannon)
HM <- Results %>% filter(Results$Cohorte == "H-M") %>% pull(Shannon)

MJ <- Results %>% filter(Results$Cohorte == "M-J") %>% pull(Shannon)
MM <- Results %>% filter(Results$Cohorte == "M-M") %>% pull(Shannon)


#Distribución normal. Aceptar H0: es dist normal
shapiro.test(Hombres)
shapiro.test(Mujeres)

shapiro.test(Jovenes)
shapiro.test(Mayores)

shapiro.test(HJ)
shapiro.test(HM)

shapiro.test(MJ)
shapiro.test(MM)

# Homogenidad
var.test(Hombres,Mujeres)
var.test(Jovenes,Mayores)
bartlett.test(Results$Shannon,Results$Cohorte) # Para dist normal. H0: varianzas son homogeneas
#levene.test(Results$Shannon,Results$Group) # Para dist no normal

boxplot(Hombres,Mujeres,names=c("Hombres","Mujeres"))
boxplot(Jovenes,Mayores,names=c("Jóvenes","Mayores"))
boxplot(HJ,HM,MJ,MM,names=c("H-J","H-M","M-J","M-M"))

#Male-Female t-test
t.test(Hombres, Mujeres, var.equal = TRUE)

#Young-Old t-test
t.test(Jovenes, Mayores, var.equal = TRUE)

#Fourties ANOVA
anova <- aov(Shannon ~ Cohorte, Results)
summary(anova)
Tukey <- TukeyHSD(anova, conf.level=0.95)
plot(Tukey, las=1 , col="brown")

#U de Mann-Whitney (Pairwise) Group
pairwise.wilcox.test(Results$Shannon, Results$Cohorte, p.adjust.method = "none" )
#--------------------------------------------------------#


#-----------RAIN-CLOUD PLOTS-----------------------------#

#After manual inspection classes with p > 0.05 in hands are:

# Aromatic Hombres jóvenesdrocarbons
# Carboximidic acids and derivatives
# Coumarins and derivatives
# Phenol esters
# Sphingolipids
# Triazines

#Aromatic hydrocarbons
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Aromatic hydrocarbons`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Aromatic hydrocarbons"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Hidrocarburos aromáticos") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Hidrocarburos aromaticos.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")


#Carboximidic acids and derivatives
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Carboximidic acids and derivatives`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Carboximidic acids and derivatives"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .2
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Ácidos carboximídicos y derivados") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Acidos carboximidicos y derivados.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")


#Coumarins and derivatives
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Dibenzylbutane lignans`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Coumarins and derivatives"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Cumarinas y derivados") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Cumarinas y derivados.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")


#Phenol esters
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Phenol esters`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Phenol esters"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Ésteres de fenol") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Ésteres de fenol.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")


#Sphingolipids
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Sphingolipids`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Sphingolipids"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Esfingolípidos") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Esfingolípidos.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")

#Triazines
p1_subc <- ggplot(T_met_statistics_Cohort, aes(x = Cohorte, 
                                               y = `Triazines`, 
                                               fill = Cohorte, colour = Cohorte,
                                               "Triazines"))  + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA,
    colour = c("#e5aecd", "#8fd9c8"),alpha = 0.4
  ) +
  scale_fill_manual(values=c("#e5aecd", "#8fd9c8")) + 
  geom_point(
    size = 5,
    aes(fill=factor(Cohorte)),
    alpha = .5,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values=c("#e5aecd", "#8fd9c8"))+
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  stat_compare_means(size = 8, label.x = 1.5, method = "t.test", 
                     label.y = c(10,19),
                     comparisons = list("Jovenes","Mayores")) +
  ggtitle("Triazinas") +
  theme(plot.title = element_text(size = 30, 
                                  hjust = 0.5, vjust = 0.5, face = "bold"))

p1_subc +
  xlab("Cohorte") + ylab("Abundancia") +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, 
                                   size = rel(1))) +
  theme(axis.text.y = element_text(size = rel(2.5))) +
  theme(axis.title = element_text(size = 4)) + #Axis Labels size
  theme(legend.text = element_text(size = 2)) + #Legend content
  theme(legend.title = element_text(size=2, face="bold")) #Legend Tittle

ggsave(".//2. Plots//Español//Triazinas.tiff",
       device = "tiff",
       width = 3000,
       height = 4025,
       units = "px",
       dpi = "retina")


