df <- read.csv("~/Downloads/Justin.csv")

library(plyr)
library(ggplot2)

df$Mass = as.numeric(as.character(df$Mass))
df = df[!is.na(df$Mass),]
df$Pupae.Type[df$Pupae.Type==" MP"]="MP"
d <- ddply(df,c("Pupae.Type","Colony_ID"),summarise,
           N = length(Mass),
           mean = mean(Mass),
           sd = sd(Mass))

d = d[d$N >= 10,]
d$CV = d$sd/d$mean

d2 <- ddply(d,c("Pupae.Type"),summarise,
            N = length(mean),
            Mean =mean(mean),
            sd = sd(mean))
d3 <- ddply(d,c("Pupae.Type"),summarise,
            N = length(CV),
            Mean =mean(CV),
            sd = sd(CV))

cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#000000", "#E69F00", "#D55E00", "#CC79A7")

ggplot(d[d$sd<0.05,],aes(x=CV,fill=Pupae.Type))+
    geom_histogram()+
    scale_fill_manual(values=cbbPalette)+
    theme_bw()

ggplot(df[df$Mass < 1,],aes(x=log(Mass),fill=Pupae.Type))+
    geom_histogram()+
    scale_fill_manual(values=cbbPalette,name="Module")+
    theme_bw()

d = read.csv("~/Downloads/SexualPupaeSamples - Sheet1.csv")
d2 <- ddply(d,"Genotype",summarise,
            Gynes=sum(Gynes),
            Males = sum(Males))
df <- data.frame(Larva = rep("JH2",6),NH = rep("JH1",6),NG=rep("JH4",6))
write.table(df,file="MockGenes.txt",row.names=FALSE)

connections <- data.frame(Type = rep(c(0,1,2),each=5),From = rep(c(3,2,6,2,4),each=3),To = rep(c(2,3,1,4,2),each=3))
write.table(connections,file="MockConnections.txt",row.names=FALSE)
