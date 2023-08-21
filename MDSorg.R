library(MVar.pt)
library(ggvegan)

mds <- planilhaorg[,1:38]
mds[,2:38] <- mds[,2:38]/orgseg[,9]
mds[,2:38] <- log(mds[,2:38]+1, base=10)
mds[,38] <- apply(mds[,2:37],1,sum) 
names(mds)[38] <- ("TotalLog")

#definir os grupos
orgseg[c(6:15,24:25),11]<- "Grupo 1"
orgseg[1:5,11]<- "Grupo 2"
orgseg[c(16:23,27:29,32:34),11]<- "Grupo 3"
orgseg[c(26,30,31,35:43),11]<- "Grupo 4"

#tabela para correlação
tabcor <- orgseg
tabcor[,12] <- log(tabcor[,10]+1, base = 10)
tabcor[,13:15] <- log(tabcor[,c(5,6,8)]+1, base =10)
tabcor[,16]<-diversity(mds[,2:38])
tabcor<-tabcor[-c(1,42,43),]


# correlação Organismos x Temperatura
cor.test(tabcor$V16,tabcor$Temperatura)
# correlação Organismos x Salinidade
cor.test(tabcor$V16,tabcor$Salinidade.1)
# correlação Organismos x Clorofila
cor.test(tabcor$V16,tabcor$Clorofila.1)

cor.test(tabcor$V16,tabcor$Profundidade)

#data <- mds[,2:37]

#Dissimilaridade segmentos
mds2 <- read.xlsx("~/TCC/ODV & R/tabelaorganismos2.xlsx",7)
data2<- t(mds2[,3:42])/1.45
colnames(data2)<-mds2[,1]

data <- mds[2:41,2:37]
row.names(data2)<-orgseg[2:41,1]
distcluster <- vegdist(data2,method = "bray")*100
grupo <- hclust(distcluster, method = "aver")
plot(grupo, hang = -1)


#dissimilaridade organismos dominantes
tabcluster <- read.xlsx("~/TCC/ODV & R/tabelaorganismos2.xlsx",4)
row.names(tabcluster)<-orgseg[,1]
tabcluster<-t(tabcluster)
tabcluster <- log(tabcluster+1)
distcluster <- vegdist(tabcluster,method = "bray")*100
grupo2 <- hclust(distcluster, method = "aver")
plot(grupo2, hang = -1)





#NMDS

dadosambiente<- tabcor[1:40,11:15]
row.names(dadosambiente)<-orgseg[1:40,1]
data2<-log(data2+1)
nmds1<- metaMDS(data2, autotransfor=FALSE)
autoplot(nmds1)
fort<-fortify(nmds1)
linhas <- envfit(nmds1$points,tabcor[,5:8],permutations = 1000 )
linhas<-fortify(linhas)

ggplot()+
  geom_point(data=subset(fort, Score=="sites"), 
             mapping = aes(x=NMDS1,y=NMDS2),
             size=4,
             alpha=0)+
  geom_segment(data=subset(fort, Score=="species"),
               mapping=aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
               arrow = arrow(length = unit(0.015,"npc"),type = "closed"),
               colour="darkgray",size=0.8)+
  geom_text(data=subset(fort, Score=="species"),
            mapping=aes(label=Label,x=NMDS1*1.5,y=NMDS2*1.5))+
  xlim(-1,1)+
  theme_bw()

linhas$vectors$arrows


ggplot()+
  geom_point(data=subset(fort, Score=="sites"), 
                   mapping = aes(x=NMDS1,y=NMDS2,color=dadosambiente$V11,shape=dadosambiente$V11),
                   size=4,
                   alpha=1)+
  geom_segment(data=linhas,
                     mapping=aes(x=0,y=0,xend=MDS1,yend=MDS2),
                     arrow = arrow(length = unit(0.015,"npc"),type = "closed"),
                     colour="darkgray",size=0.8)+
  geom_text(data=linhas,
                  mapping=aes(label=Label,x=MDS1*1.2,y=MDS2*1.2))+
  geom_text(data=subset(fort, Score=="sites"),
            mapping=aes(label=Label,x=NMDS1*1.1,y=NMDS2*1.1))+
  xlim(-1,1)+
  theme_bw()




#
md <- MDS(data = distcluster, distance = "euclidean", title = NA, xlabel = NA,  
          ylabel = NA, posleg = 2, boxleg = TRUE, axes = TRUE, color = TRUE,  
          linlab = NA, class = cls,
          savptc = FALSE, width = 3236, height = 2000, res = 300)

print("Matriz das distancias:"); md$mtxD

grupo <- hclust(dist(data))


#Diversidade e equitabilidade

ggplot(mds, aes(x=Segmento,y=diversity(mds[,2:37]))) + geom_point(color="darkblue") +geom_line(color="darkblue") +
  scale_x_continuous(breaks = seq(from = 1, to = 87, by = 2)) +
  geom_line(aes(mds[,1],((diversity(mds[,2:37]))/(log(length(which(mds[,2:37]!=0))))*6)),color="red" ) +
  geom_point(shape = 15, color= "red",aes(mds[,1],((diversity(mds[,2:37]))/(log(length(which(mds[,2:37]!=0))))*6)))+
  scale_y_continuous(sec.axis = sec_axis(~./6, name = "Equitabilidade")) +
  labs(x="Segmento CPR", y="Diversidade (bits.ind-1)") +theme_bw()

