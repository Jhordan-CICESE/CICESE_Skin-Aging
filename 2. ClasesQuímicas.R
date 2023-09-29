library(dplyr)
library(plotly)
#Plotly
Sys.setenv("plotly_username" = "froz9")
Sys.setenv("plotly_api_key" = "3YcsdCmBtC8rfYLZA3PR")

#SET ENVIRONMENT
setwd("G://Mi unidad//CICESE//Tesis//2. Analysis//2. Results//2. FBMN//2. ChemicalClasses")

#Get Class Colors
colors <- read.csv(paste(".//0. Source////ColorCode.csv", sep = ""))
colors <- colors[,-1] # -1 = Español, -2 = English
colnames(colors)[1] <- "Class"

# Set Name Variable
name <- "FBMN"
description <- "Únicos y filtrados por error de masa < 10 ppm"
# Sin filtrar / FBMN
# Filtrados por error de masa < 10 ppm / _Filtered
# Únicos y filtrados por error de masa < 10 ppm / Filtered_Unique
# Grupo de jóvenes / _Young
# Grupo de mayores / _Old
filter <- 6

#Read Table
Canopus <- read.csv(paste(".//0. Source//Classes_",name,".csv", sep = ""))
class <- data.frame(table(Canopus$Clase)) %>% rename(Clase = Var1)  
# Clase = Español, Class = English
colnames(class)[1] <- "Class"

df1 <- class %>% filter(Freq < filter)
df2 <- class %>% filter(Freq >= filter)
df3 <- as.data.frame(merge("Otros con < 3",sum(df1$Freq)))  %>% 
  rename(Class = x, Freq = y)
df4  <- rbind(df2,df3)
df4 <- merge(colors,df4, by = "Class") 
X <- df4[df4$Class %in% c("Otros con < 3","Sin asignar"),]
df4 <- df4[df4$Class != "Sin asignar", ]
df4 <- df4[df4$Class != "Otros con < 3", ]
df4 <- rbind(df4,X)

colnames(df4) <- c("Class","Color", "Freq")

fig_class <- plot_ly(df4, labels = ~Class, values = ~Freq, type = 'pie', 
                     sort = FALSE, 
                     direction = "clockwise",
                     marker = list(colors= ~Color), 
                     textfont = list(size = 10))%>%
  layout(font=list(size = 10),
         legend = list(title=list(text=paste('<b> Clases químicas de Canopus\n',description)), 
                       x= 100, y=0.5))
fig_class

plotly_IMAGE(fig_class, format="png", 
             width = 1000, height = 1000, 
             out_file=paste(".//1.2 Resultados//Class_",
                            name,"_",filter,".png", sep = ""))
