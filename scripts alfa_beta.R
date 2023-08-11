
install.packages("iNEXT")
install.packages("ggplot2")
install.packages("RColorBrewer")

library(iNEXT)
library(ggplot2)
library(RColorBrewer)

getwd()
setwd("C:/Users/camil/OneDrive - unillanos.edu.co/2022-II/Taxonomía animal/Módulo de peces/Proyecto/Proyecto_R_div_alfa_beta")

### Diversidad alfa

divalfa <- read.table("divalfa.csv", header = T, sep = ";", row.names = 1)
divalfa

View(divalfa)

diversidadq012 <- iNEXT(divalfa, q=c(0,1,2), datatype="abundance", size=NULL,
                      endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=100)
diversidadq012$DataInfo
View(diversidadq012$DataInfo)

#q0

diversidadq0 <- iNEXT(divalfa, q=0, datatype="abundance", size=NULL,
                      endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=100)


ggiNEXT(diversidadq0, type = 3, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) + ylim(c(0,23))+   xlim(c(0.25, 1))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=0.75, y=16, label= "q=0", cex=7) +
  labs(x = "Cobertura del muestreo", y = "Riqueza (q=0)")

ggiNEXT(diversidadq0, type = 1, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=100, y=24, label= "q=0", cex=7) +
  labs(x = "Número de individuos", y = "Riqueza (q=0)")

#q1

diversidadq1 <- iNEXT(divalfa, q=1, datatype="abundance", size=NULL,
                      endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=100)

ggiNEXT(diversidadq1, type = 3, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) + ylim(c(0.5,6))+   xlim(c(0.25, 1))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=0.4, y=4, label= "q=1", cex=7) +
  labs(x = "Cobertura del muestreo", y = "Shannon-Wienner (q=1)")

ggiNEXT(diversidadq1, type = 1, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) + theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=300, y=3, label= "q=1", cex=7) +
  labs(x = "Número de individuos", y = "Shannon-Wienner (q=1)")

#q2

diversidadq2 <- iNEXT(divalfa, q=2, datatype="abundance", size=NULL,
                      endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=100)

ggiNEXT(diversidadq2, type = 3, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) + ylim(c(0.5,4))+   xlim(c(0.25, 1))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=0.4, y=3, label= "q=2", cex=7) +
  labs(x = "Cobertura del muestreo", y = "Dominancia de Simpson (q=2)")

ggiNEXT(diversidadq1, type = 1, facet.var= "Order.q", color.var = "Order.q") + 
  theme_bw(base_size = 8) + theme(legend.position = "bottom",legend.title = element_blank())+
  annotate("text", x=300, y=2, label= "q=2", cex=7) +
  labs(x = "Número de especies", y = "Dominancia de Simpson (q=2)")

### Diversidad beta

install.packages("betapart")

library(betapart)
library(ggplot2)

divbeta <- read.table("pres_aus.csv", header = T, sep = ";", row.names = 1)
divbeta
View(divbeta)

ceram.s.core <- betapart.core(divbeta)
ceram.s.core

ceram.s.multi <- beta.multi(ceram.s.core)
ceram.s.multi
#=
abundancia.multi <- beta.multi.abund(divbeta,index.family="bray")
abundancia.multi

ceram.s.samp <- beta.sample(ceram.s.core, sites=5, samples=100)
ceram.s.samp
dist.s <- ceram.s.samp$sampled.values

plot(density(dist.s$beta.SOR),xlim=c(0,1), ylim=c(0,20), xlab="Beta diversity", main="", lwd=3)
lines(density(dist.s$beta.SNE), lty=1, lwd=2) 
lines(density(dist.s$beta.SIM), lty=2, lwd=2)

pair.s <- beta.pair.abund(divbeta)
pair.s
pair.s$beta.bray.bal
pair.s$beta.bray.gra
pair.s$beta.bray

#plot

dist.s <- ceram.s.samp$sampled.values
plot(hclust(pair.s$beta.bray.bal, method="average"), hang=-1, main="", sub="", xlab="")
title(xlab=expression(beta[bray.bal]), line=0.4)

plot(hclust(pair.s$beta.bray.gra, method="average"), hang=-1, main="", sub="", xlab="")
title(xlab=expression(beta[bray.gra]), line=0.4)

# Beta Presencia-Ausencia

install.packages("permute")
install.packages("lattice")
install.packages("vegan")

library(permute)
library(lattice)
library(vegan)
 
betapa <- read.table("pres_aus.csv", header = T, sep = ";", row.names = 1)
betapa
View(betapa)

# decostand es para estandarizar y dejar todo de 0 y 1

betaps <- decostand(betapa, "pa")
betaps

# Correrlos con bata core

presencia <- betapart.core(betaps)

presencia.multi <- beta.multi(betaps,index.family="jaccard")
presencia.multi

presencia.samp <- beta.sample(presencia,index.family="jaccard", sites=5, samples=100)
presencia.samp 

dist <- ceram.s.samp$sampled.values
dist

pair.presencia<- beta.pair(betaps)
pair.presencia



dist. <- ceram.s.samp$sampled.values
plot(hclust(pair.presencia$beta.sim, method="average"),hang=-1, main='', sub='', xlab='', ylab= "Disimilitud")
title(xlab=expression(beta[sim]), line=0.7)
plot(hclust(pair.presencia$beta.sne, method="average"),hang=-1, main='', sub='', xlab='', ylab= "Altura")
title(xlab=expression(beta[sne]), line=0.7)
