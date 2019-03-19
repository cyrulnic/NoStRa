library(ggplot2)

# DATA
#set.seed(345)
#Sector <- rep(c("S01","S02","S03","S04","S05","S06","S07"),times=7)
#Year <- as.numeric(rep(c("1950","1960","1970","1980","1990","2000","2010"),each=7))
#Value <- runif(49, 10, 100)
#data <- data.frame(Sector,Year,Value)
#ggplot(data, aes(x=Year, y=Value, fill=Sector)) +   geom_area()


#Sector <- rep(c("S01","S02","S03"),times=4)
#Year <- as.numeric(rep(c("1950","1960","1970","1980"),each=3))
#Value <- runif(12, 10, 100)
#data <- data.frame(Sector,Year,Value)
#ggplot(data, aes(x=Year, y=Value, fill=Sector)) +   geom_area()

# data go top down abscissa-category by abscissa-category:
# Zero: Be, Br, Bc
# Weak: Be, Br, Bc
# ...


Sector = rep(c("w_e","w_es","w_r","w_c"),times=4)
Nonstationarity_Strength = rep(c("Zero","Weak","Medium","Strong"),each=4)
Nonstationarity_Strength = as.numeric(rep(c("0","1","2","3"),each=4))
#Weights = c(0,0,1,  0.2,0.64,0.16,  0.2,0.72,0.08,  0.2,0.8,0)

# weights of the opt HHBEF: 4 columns = 4 nstatio strengths

w_cvr = c(1, 0.07, 0.09,  0)
w_evp = c(0, 0.1,  0.20,  0.2)
kappa0= c(1, 0.25, 0.25, 1)

w=1-w_cvr
mu=1-w_evp

# $w_{\rm e}  = (1-\mu)\kappa_0$ 
# $w_{\rm es} = (1-\mu)(1-\kappa_0)$ 
# $w_{\rm c}  = {\mu (1- w)}/{(1 -\mu w)}$ 
# $w_{\rm r}  = \mu w(1-\mu)/(1-\mu w)$ 

w_e  = (1-mu)*kappa0
w_es = (1-mu)*(1-kappa0)
w_r  = (1-mu) * mu*w/(1 - (mu*w))
w_c  = mu * (1-w)/(1 - (mu*w))

# sanity check (sh.be all ones)
w_e+w_es+w_c+w_r

Weights = rbind(w_e, w_es, w_r, w_c)
Weights
apply(Weights, 2, sum) # sh.be all ones
weights=as.vector(Weights)

data = data.frame(Sector, Nonstationarity_Strength, weights)
p=ggplot(data, aes(x=Nonstationarity_Strength, y=weights, fill=Sector),
         cex.main=1.7, cex.axis=1.3, cex.lab=1.6) +  
  geom_area()
#geom_area(position = position_stack(reverse = T))

#pngname="CVM_weights.png"
#png(pngname, width=6, height=6, units = "in", res=300)
#par(mai=c(1.2,1.2,0.7,0.7))
#
namefile=paste0("CVM_weights.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
  p + ggtitle("Weights of B_c, B_e, B_es, B_r") +
    theme(text = element_text(size=16))
dev.off()
           

#p

#p + scale_fill_grey(start = 0.7, end = 0.3)

namefile=paste0("CVM_weights_brewer.pdf")
pdf(namefile, width=7, height=7)
par(mgp=c(2.5, 1, 0))
p + ggtitle("Weights of B_c, B_e, B_es, B_r")  +
  theme(text = element_text(size=16)) + 
    scale_fill_brewer(palette = 1)
dev.off()
            
#interchange order : reverse

#Sector = rep(c("B_clim","B_recent","B_ensm"),times=4)
#Nonstationarity_Strength = rep(c("Zero","Weak","Medium","Strong"),each=3)
#Nonstationarity_Strength = as.numeric(rep(c("0","1","2","3"),each=3))
#Value = c(1,0,0,  0.16,0.64,0.2,  0.08,0.72,0.2,  0,0.8,0.2)
#data = data.frame(Sector, Nonstationarity_Strength, Value)
#ggplot(data, aes(x=Nonstationarity_Strength, y=Value, fill=Sector)) +  
#  geom_area(position = "stack")
#  #geom_area(position = position_stack(reverse = T))

         
