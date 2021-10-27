library(ggplot2)
library(tidyverse)
######################
# HEAT MAP
######################
#---- for full complexity ----
data.final <- NULL
for(FocalPrefix in c("LEMA","CHFU","HOMA","MESU")){
  data.final.focal <- NULL
data.alpha <- read.csv(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_ij_",FocalPrefix,".csv"))
data.alpha<- data.alpha[,-1]

if (FocalPrefix == "HOMA" | FocalPrefix == "MESU" | FocalPrefix == "CHFU" ){
  data.final.focal  <- data.frame(year=data.alpha[,ncol(data.alpha)])
  data.alpha <- data.alpha[,-ncol(data.alpha)]
}else{data.final.focal  <- data.frame(year=years)}

data.final.focal$inter_alpha <- rowSums(data.alpha)

data.beta <- read.csv(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/beta_Inclusion_",FocalPrefix,".csv"))
data.beta <- data.beta[,-1]
names(data.beta) <- c(get(paste0("SpNames_",FocalPrefix)),"year")
if (FocalPrefix == "LEMA"){data.beta$year <- rep(years,each=length(get(paste0("SpNames_",FocalPrefix))))}
data.beta<-  data.beta%>%
  group_by(year) %>%
  summarize_all(sum) 
data.beta <- as.data.frame(data.beta)

data.final.focal$focal <- FocalPrefix
data.final.focal$intra_beta <- data.beta[,FocalPrefix]
data.final.focal$inter_beta <- rowSums(data.beta[,which(!names(data.beta) %in% c(FocalPrefix,"year"))])

data.final.focal  <- gather(data.final.focal,inter_alpha  , inter_beta ,intra_beta,
                            key = "interaction.type", value = "number")
data.final$year <- as.numeric(data.final$year)
data.final.focal$year <- as.numeric(data.final.focal$year)
data.final  <- bind_rows(data.final,data.final.focal)

}

# addition of all possibility
expansion.final.data <- expand.grid(unique(data.final$year),
                                    unique(data.final$focal),
                                    unique(data.final$interaction.type))
names(expansion.final.data ) <- c("year","focal",'interaction.type')
str(expansion.final.data)
expansion.final.data <- as.data.frame(expansion.final.data)
expansion.final.data$year <- as.numeric(expansion.final.data$year )
expansion.final.data$focal <- as.character(expansion.final.data$focal)
expansion.final.data$interaction.type <- as.character(expansion.final.data$interaction.type)
data.final <- dplyr::left_join(expansion.final.data ,data.final,
                               by = c("year","focal",'interaction.type'))

data.final$complexity <- c("Species")
data.final.full <-data.final
plot.map <- list()
for (interactions.levels in levels(as.factor(data.final$interaction.type))) {
  if(interactions.levels == "inter_alpha"){ plot.title ="Alpha inter-"} 
  if(interactions.levels == "inter_beta"){ plot.title ="Beta inter-"} 
  if(interactions.levels == "intra_beta"){  plot.title ="Beta intra-"} 
  
  
  plot.map[[paste("plot",interactions.levels , sep = "_")]] <- ggplot(data = data.final[which(data.final$interaction.type == interactions.levels),],
                                                                      aes(x=year, y=focal, fill=number)) + 
    geom_tile() +
    geom_text(aes(label=number),color="white") + 
    ggtitle(plot.title)+
    scale_fill_gradientn(colours = c("black", "blue", "red","yellow"), 
                         #values = c(0,0.01,1),
                         limits= c(0,20)) +
    theme_black_noline() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())

}
ggarrange(plotlist = plot.map,nrow = 1, ncol=3, 
          common.legend = T,legend= "bottom")  
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/heatplot_full",".pdf"),
  plot = last_plot())

#---- for less complexity ----
data.final <- NULL
for (FocalPrefix in c("MESU","LEMA","HOMA","CHFU")){
  for(complexity  in c("family","class","rareORabundant")){

data.alpha<- read.csv(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_ij_",complexity,"_",FocalPrefix,".csv"))
data.alpha<- data.alpha[,-1]
data.alpha.final <- NULL

data.alpha.final$year <-data.alpha[,"year" ]
data.alpha.final$inter_alpha <- rowSums(data.alpha[,which(!names(data.alpha) %in% c(FocalPrefix,"year"))])

data.alpha.final <- as.data.frame(data.alpha.final)
data.alpha.final$focal <- FocalPrefix
data.alpha.final$complexity <- complexity

    
data.beta<- read.csv(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/beta_Inclusion_family_",FocalPrefix,".csv"))
data.beta<- data.beta[,-1]
data.beta<-  data.beta%>%
   group_by(year) %>%
   summarize_all(sum) 
data.beta <- as.data.frame(data.beta)

data.beta.final <- data.frame(year=data.alpha.final$year)
data.beta.final$focal <- FocalPrefix
data.beta.final$complexity <- complexity

data.beta.final$intra_beta <- data.beta[,FocalPrefix]
data.beta.final$inter_beta <- rowSums(data.beta[,which(!names(data.beta) %in% c(FocalPrefix,"year"))])

data.final.focal  <- gather(right_join(data.alpha.final,data.beta.final, 
                                 by=c("year","focal","complexity")),
                      inter_alpha  , inter_beta ,intra_beta,
       key = "interaction.type", value = "number")

data.final  <- bind_rows(data.final,data.final.focal)

}
}
#data.final$year <- as.factor(data.final$year)

# addition of all possibility
expansion.final.data <- expand.grid(unique(data.final$year),
                                     unique(data.final$focal),
                                     unique(data.final$complexity),
                                     unique(data.final$interaction.type))
names(expansion.final.data ) <- c("year","focal","complexity",'interaction.type')

data.final <- dplyr::left_join(expansion.final.data ,data.final,
                                 by = c("year","focal","complexity",'interaction.type'))

# make the  plot
for (complexity in c("family","class","rareORabundant")){
  plot.map <- list()
for (interactions.levels in levels(as.factor(data.final$interaction.type))) {
  if(interactions.levels == "inter_alpha"){ plot.title ="Alpha inter-"} 
  if(interactions.levels == "inter_beta"){ plot.title ="Beta inter-"} 
  if(interactions.levels == "intra_beta"){  plot.title ="Beta intra-"} 

plot.map[[paste("plot",interactions.levels , sep = "_")]] <- ggplot(data = data.final[which(data.final$interaction.type == interactions.levels &
                                                                                              data.final$complexity==  complexity),],
                                                                    aes(x=year, y=focal, fill=number)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("black", "blue", "red","yellow"), 
                      # values = c(0,0.01,1),
                       limits= c(0,20)) +
  geom_text(aes(label=number),color="white") + 
  theme_black_noline() +
  ggtitle(plot.title) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
}
#pdf(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/heatplot_",complexity,".pdf"))
ggarrange(plotlist = plot.map,nrow = 1, 
          ncol=3, common.legend = T, legend="bottom")  
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/heatplot_",complexity,".pdf"),
  plot = last_plot())
#dev.off()
}
library(ggpubr)

######################
# LINE graphs 
######################
# ---- number of interactions ----
str(data.final)

data.final$year <- as.numeric(data.final$year)
data.final.full$year <- as.numeric(data.final.full$year)
data.final.all <- bind_rows(data.final.full,data.final)
data.final.all$complexity <- factor(as.factor(data.final.all$complexity), 
                                          levels = c("rareORabundant", "class", "family","Species"))

cbPalette <- cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 ggplot(aggregate(number ~ interaction.type + complexity,
                                    g.interaction.type, sum), 
       aes(x=complexity , y= number/3, 
           #color=interaction.type,
           fill=interaction.type,
           group=interaction.type)) + 
   geom_bar(stat='identity',position="fill") + 
   scale_fill_manual(values=cbPalette) +
  #geom_line(size=2) +
  #geom_point(size=3) + 
  scale_colour_manual(values=cbPalette) +
  #theme_black() + 
   theme_bw() + 
  ylab("number of interactions") + 
  scale_x_discrete(labels=c("abundance","forb or grass",
                          "family","species")) +
  
  theme(axis.title.x = element_blank()) +
  guides(color=guide_legend(title="")) 
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_complexity_interaction.pdf"),
  plot = last_plot())

 ggplot(aggregate(number ~ complexity + year + focal,
                                    g.interaction.type, sum), 
                          aes(x=focal , y= number/3, 
                              #color=complexity,
                              fill=complexity,
                              group=complexity)) + 
   geom_bar(stat='identity',position="fill") + 
   scale_fill_manual(values=rev(cbPalette)) +
  #geom_line(size=2) +
  #geom_point(size=3) + 
  scale_colour_manual(values=rev(cbPalette)) +
  #theme_black() + 
   theme_bw() + 
  ylab("number of interactions") + 
   scale_x_discrete(labels=c("abundance","forb or grass",
                             "family","species")) +

  theme(axis.title.x = element_blank()) +
  guides(color=guide_legend(title="")) 
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_complexity_interaction.pdf"),
  plot = last_plot())


year.plot <- ggplot(aggregate(number ~interaction.type + year, 
                              data.final.all, sum), 
       aes(x=year , y= number,
           fill=interaction.type,
           group=interaction.type))+ 
  geom_bar(stat='identity',position="fill") + 
  scale_fill_manual(values=cbPalette) +
  #geom_line(size=2) +
  #geom_point(size=3) + 
  scale_colour_manual(values=cbPalette) +
 #theme_black() + 
  theme_bw()+
  ylab("number of interactions") + 
  scale_x_continuous(breaks=c(2017,2018,2019)) +
  theme(axis.title.x = element_blank()) +
  guides(color=guide_legend(title="")) 

ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_year_interaction.pdf"),
  plot = last_plot())



species.plot.complexity <- ggplot(aggregate(number ~ focal + complexity,
                                            data.final.all, sum), 
                    aes(x=focal , y= number, 
                        fill=complexity,
                        group=complexity)) + 
  geom_bar(stat='identity',position="fill") + 
  scale_fill_manual(values=rev(cbPalette)) +
  #geom_line(size=2) +
  #geom_point(size=3) + 
  scale_colour_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme_bw() + 
  ylab("number of interactions") + 
  #scale_x_discrete(labels=c("2017","2018","2019")) +
  theme(axis.title.x = element_blank()) +
  guides(color=guide_legend(title="")) 
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_complexity_species.pdf"),
  plot = last_plot())

ggarrange(complexity.plot ,year.plot,
          common.legend = T,legend="bottom")
ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/line_year.pdf"),
  plot = last_plot())


# ---- percentage of interactions ----
    g.interaction.type <- NULL
    n.complexity <- NULL
for (n in levels(as.factor(data.final.all$complexity))  ) {
  g <- data.final.all[which(data.final.all$complexity ==n),]
  if (n =="Species")( n <- tolower(n))
  l <- nlevels(as.factor(plant.class[,n]))
  for (i in levels(as.factor(data.final.all$interaction.type))){
    k <- g[which(g$interaction.type==i),]
    if ( i == "inter_alpha"){
    if (n =="Species"){
      k$interaction.type <- 'Alpha inter'
        for (sp in levels(as.factor(data.final.all$focal))){
        r <- k[which(k$focal==sp),]
        m <- length(get(paste0("SpNames_",sp)))}}else{
          k$interaction.type <- 'Alpha inter'
          m <- l 
        }
    }
    if ( i == "inter_beta"){
      if (n =="Species"){
        k$interaction.type <- 'HOIs inter'
        for (sp in levels(as.factor(data.final.all$focal))){
        r <- k[which(k$focal==sp),]
        m <- length(get(paste0("SpNames_",sp)))^2 }}else{
          k$interaction.type <- 'HOIs inter'
          m <- l*l - 1 
        }}
    if ( i == "intra_beta"){     k$interaction.type <- 'HOIs intra'
    m <- 1}

    k$number <- k$number/m
    g.interaction.type  <- bind_rows(g.interaction.type,k)
  }
}
view(  g.interaction.type)
levels(as.factor(g.interaction.type$complexity.2))
g.interaction.type <- g.interaction.type %>%
        mutate(complexity.2=case_when(complexity=="rareORabundant" ~ "rare or abundant",
                               complexity=="class" ~ "forb or grass",
                               complexity=="family" ~ "family level",
                               complexity=="Species" ~ "species level"))
g.interaction.type$complexity.2 <-   factor(as.factor(g.interaction.type$complexity.2), 
       levels = c("rare or abundant", "forb or grass",
                  "family level","species level"))

# bar :  year interaction
ggplot(aggregate(number ~ year + interaction.type ,
                 g.interaction.type,sum), 
           aes(x=year , y= (number/4)/4, 
               fill=  interaction.type ,
               #color=  interaction.type ,
               #shape=as.factor(complexity),
               group =    interaction.type)) + 
  geom_bar(stat='identity',position="fill") + 
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(name="Interaction",
                     values=rev(cbPalette)) +
  geom_line(size=1.5,color="black") +
  geom_point(size=4,aes(shape=interaction.type )) + 
      #scale_colour_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      #scale_colour_manual(values=cbPalette) +
      theme_bw() + 
      ylab("percentage of interactions") + 
      theme(axis.title.x = element_blank(),
            legend.title = element_text(color = "white")  
      ) +
      guides(color=guide_legend(title="")) 

ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_year_interaction.pdf"),
  plot = last_plot())

# bar :  focal complexity
ggplot(aggregate(number ~ focal + complexity.2 ,
                 g.interaction.type,sum), 
       aes(x=focal , y= (number/3)/3, 
           fill= complexity.2 ,
           #shape=as.factor(complexity.2),
           group =   complexity.2)) + 
  geom_bar(stat='identity',position="fill") + 
  scale_fill_manual(values=rev(cbPalette)) +
  geom_line(size=1.5) +
  geom_point(size=4,aes(shape= complexity.2)) + 
  scale_colour_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  
  #scale_colour_manual(values=cbPalette) +
  theme_bw() + 
  ylab("percentage of interactions") + 
  
  #scale_x_discrete(labels=c("abundance","forb or grass",
  #                          "family","species")) +
  
  theme(axis.title.x = element_blank(),
        legend.title = element_text(color = "white")) +
  guides(color=guide_legend(title="")) 

ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_sp_complexity.pdf"),
  plot = last_plot())

ggplot(aggregate(number ~ focal + interaction.type ,
                 g.interaction.type,sum), 
       aes(x=focal , y= (number/3)/4, 
           fill=  interaction.type ,
           #shape=as.factor(complexity.2),
           group =    interaction.type)) + 
  geom_bar(stat='identity',position="fill") + 
  geom_line(size=1.5) +
  geom_point(size=4,aes(shape=interaction.type)) + 
  #scale_colour_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  
  scale_fill_manual(values=cbPalette) +
  theme_bw() + 
  ylab("percentage of interactions") + 
  
  #scale_x_discrete(labels=c("abundance","forb or grass",
  #                          "family","species")) +
  
  theme(axis.title.x = element_blank(),
        legend.title = element_text(color = "white")) +
  guides(color=guide_legend(title="")) 

ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_sp_interaction.pdf"),
  plot = last_plot())


ggplot(aggregate(number ~ interaction.type + complexity.2,
                 g.interaction.type,sum), 
       aes(x=complexity.2 , y= (number/4)/3, 
           fill=  interaction.type ,
           #shape=as.factor(complexity.2),
           group =    interaction.type)) + 
  geom_bar(stat='identity',position="fill") + 
  scale_fill_manual(values=cbPalette) +
  geom_line(size=1.5) +
  geom_point(size=4, aes(shape=interaction.type)) + 
  #scale_colour_manual(values=c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  
  scale_colour_manual(values=cbPalette) +
  theme_bw() + 
  ylab("percentage of interactions") + 
  
  #scale_x_discrete(labels=c("abundance","forb or grass",
  #                          "family","species")) +
  
  theme(axis.title.x = element_blank(),
        legend.title = element_text(color = "white")) +
  guides(color=guide_legend(title="")) 

ggsave(
  paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/bar_complexity_interaction.pdf"),
  plot = last_plot())
