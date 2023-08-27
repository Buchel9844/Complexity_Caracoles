gens = 1
state = c(0,1)
  params <- list(
    # stable parameters
    sim= 1, 
    g = c(gi=1,gj=1), # germination rate of seeds
    s = c(si= 1,sj= 1), # survival rate of seeds
    lambda = c(lambdai= 1 ,lambdaj= 1), # intrinsic fecundity
    # variable parameters
    
    a = matrix(ncol=2, nrow = 2, # a_slope, only positive values
                       data = c(0.1,0.1,0.1,0.1),
                       dimnames = list(c("i", "j"),
                                       c("i", "j"))),
  )

# evaluate the magnitude of interaction on fecundity
Ricker_solution <- function(gens,
                            state,
                            pars) {
  g <- pars$g # germination rate 
  s <- pars$s #seed survival
  lambda <- pars$lambda # intrinsic growth rate
  a_initial <- pars$a_initial # which int.function

  df <- data.frame( t=0:gens,  Ni=numeric(1+gens),  Nj =numeric(1+gens) )
  df[1,2:3] <- c(state[1],state[2]) #species i initial densities
  
  for(t in 1:gens){
    
    Ni <- df[t,"Ni"] # species i densities
    Nj <- df[t,"Nj"] # species j  densities

      aii <- a_initial[1,1]
      aij <- a_initial[1,2]
      aji <- a_initial[2,1]
      ajj <- a_initial[2,2]

    Fi <-  exp(lambda[1] + aii * g[1]*Ni + aij *g[2]*Nj)
    Fj <-  exp(lambda[2] + ajj * g[2]*Nj + aji *g[1]*Ni)
    
    Nit1 <- ((1-g[1]) * s[1] + g[1] * Fi)*Ni
    Njt1 <- ((1-g[2]) * s[2] + g[2] * Fj)*Nj
    df[t+1,2:3] <- c(Nit1, Njt1)
  }
  names(df) <- c("time","Ni","Nj")
  return(df)
}


test <- data.frame(x=seq(1,length(seq(-1,1,0.01)),1),
                   alpha=seq(-1,1,0.01))
lambda  = 4
test <- test %>%
  mutate(fec =exp(lambda + test$alpha ),
         col.val = case_when((- 0.1 < alpha  &  alpha  <  0.1) ~ 1,
                              (- 0.25 < alpha  &  alpha <= - 0.1) | (0.1 <= alpha  &  alpha < 0.25) ~ 2,
                              (- 0.5 < alpha  &  alpha <= - 0.25 )| (0.25 <= alpha &  alpha  < 0.5) ~ 3,
                              (alpha <= - 0.5 )| (alpha >=  0.5) ~ 4))
 
head(test)

Oneunit <- (max(test$fec) - min(test$fec))/100

fec = function(x){
  exp(lambda + x)
}

fec.proportion = function(x){
  (x-exp(lambda))/(Oneunit*100)
}
fec.proportion.rev = function(x){
  fec =exp(lambda  + x)
  fec.prop = (fec -exp(lambda))/(Oneunit*100)
  return(fec.prop)
}


#cumsum((fec))
breaks.vec <- c(-1,-0.75,-0.5,-0.25,-0.1,0,0.1,0.25,0.5,0.75,1)
breaks.y.vec <- fec(breaks.vec)
breaks.y.lab <- fec(breaks.vec)

plot.scale.fec <- ggplot(test,aes(x=alpha, y=fec, fill= as.factor(col.val))) +   
  geom_col(na.rm=TRUE) +
  scale_y_continuous(
    "Fecundity",
    breaks=fec(breaks.vec),
    labels=round(fec(breaks.vec)),
    sec.axis = sec_axis(fec.proportion, name = "Proportion of fecundity change",
                        labels=paste0(round(fec.proportion.rev(breaks.vec)*100),"%"),
                        breaks=fec.proportion.rev(breaks.vec))
  ) +
  scale_fill_manual("",values= c("#999999" ,"#56B4E9","#F0E442","#CC79A7"),
                         labels=c("insignificant","small","medium","large"))+
  theme_bw() +
  scale_x_continuous("interaction effect", expand=c(0.01,0.01))+
  coord_cartesian(xlim = c(-1, 1), clip="off") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 50, b = 0, l = 0))) +
  annotate("text", x=-1.1 ,y=exp(lambda), size=4,
           parse = TRUE, label ="e^(lambda)") +
  annotate("text", x=-1.15 ,y=exp(lambda + 0.25), size=4,
           parse = TRUE, label ="e^(lambda+0.25)") +
  annotate("text", x=-1.15 ,y=exp(lambda - 0.25), size=4,
           parse = TRUE, label ="e^(lambda-0.25)")

for( n in breaks.vec){
  plot.scale.fec <-  plot.scale.fec +
    annotate("segment", x=n ,xend=1.03,
             y=exp(lambda + n ) ,yend=exp(lambda+ n ) ,
             linetype="dotted",
             col="black")
}


test.density <- data.frame(alpha=seq(-1,1,0.001)) %>%
  mutate(gaussian = dnorm(alpha,mean = 0, sd = 0.5, log = FALSE),
         col.val = case_when((- 0.1 < alpha  &  alpha  <  0.1) ~ 1,
                    (- 0.25 < alpha  &  alpha <= - 0.1) | (0.1 <= alpha  &  alpha < 0.25) ~ 2,
                    (- 0.5 < alpha  &  alpha <= - 0.25 )| (0.25 <= alpha &  alpha  < 0.5) ~ 3,
                    (alpha <= - 0.5 )| (alpha >=  0.5) ~ 4))


plot.density <- ggplot(test.density,aes(y=gaussian,x=alpha,fill=as.factor(col.val))) +
  geom_col(na.rm=TRUE) +
  scale_fill_manual("",values= c("#999999" ,"#56B4E9","#F0E442","#CC79A7"),
                    labels=c("insignificant","small","medium","large"))+
  scale_x_continuous("interaction effect",
                     expand=c(0.01,0.01)) +
  labs(y="")+
  theme_bw() +
  coord_cartesian(xlim = c(-1, 1), clip="off") 


ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/scale_alpha.pdf"),
       dpi="retina",scale=1, 
       ggarrange(plot.scale.fec ,plot.density ,ncol=1,nrow=2,heights = c(3,1),
                 align = c("v"),
                 common.legend =T,legend="bottom")
       )





