

Abundance <- read.csv("/Users/lisabuche/Code/Project/Caracoles/UnderTheHood/data/abundance.csv")
str(Abundance)
library(ggplot2)
AbundanceFocal <- subset(Abundance , species %in% c("LEMA","CHFU","CETE","HOMA"))
str(AbundanceFocal)
AbundanceFocal <- aggregate(individuals ~ month + species, 
                            data= Abundance,mean)

AbundanceFocal$date <- paste(AbundanceFocal$year,AbundanceFocal$month, sep="_")
ggplot(AbundanceFocal, aes(y=individuals,x=month)) +
  geom_bar(stat="identity",aes(fill=species),position = "dodge") +
  ylab("individuals per 1m square plot") + scale_y_sqrt() + 
  geom_smooth(method = "lm")
  

