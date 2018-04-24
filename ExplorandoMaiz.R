##Este script explora los datos de maizteocitle
rm(list=ls())
#Leer los datos
maices<-read.delim("../meta/maizteocintle_SNP50k_meta_extended.txt",header = T)
#Que tipo de objeto se creo
class(maices)
#Ver las primers 6 filas de los datos 
head_maices<-maices[c(1:6),]
View(head_maices)
#Cuántas muestras hay
nrow(maices)
#De cuántos estados hay muestras
nlevels(maices$Estado)
#Cuántas muestras fueron colectadas antes de 1980 
length(which(maices$A.o._de_colecta<1980))
#Cuántas muestras hay de cada raza
summary(maices$Raza)
#promedio de altitud
mean(maices$Altitud)
#altura máxima y mínima
max(maices$Altitud)
min(maices$Altitud)
#Crea nueva data frame con sólo las muestras de olotillo
olotillo<-maices[which(maices$Raza=="Olotillo"),]
#Crea nueva data frame con las muestras de Reventador, Jala y Ancho
Reventador_Jala_Ancho<-maices[which(maices$Raza=="Reventador" | maices$Raza=="Jala" | maices$Raza=="Ancho"),]
#Guardar nueva data frame
write.csv(Reventador_Jala_Ancho,"../meta/submat.csv")
