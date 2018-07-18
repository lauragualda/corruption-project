library(MASS)
df <- read.csv("n_irreg.csv", sep=";", header = T)

# estimate NB(r,p)
theta <- fitdistr(df$n_irreg_EDUC[which(df$n_irreg_EDUC>0)]-1,densfun = "negative binomial")
r = theta$estimate[["size"]]
p = theta$estimate[["size"]]/sum(theta$estimate)

# function that calculates E q^(N+1) where N~NB(r,p), 0<q<1
Eqn<-function(q, r, p){q*((1-p)/(1-p*q))^r}

# goodness of fit
O<-c()
nobs <- sum(df$n_irreg_EDUC>0)
nmax <- max(df$n_irreg_EDUC)-1
for(i in 0:nmax){O[i+1] <- sum(df$n_irreg_EDUC==i+1)/nobs}

E<-dnbinom(seq(0,nmax),r,p)
E<-E/sum(E) #normalize 

# q-q plot
plot(x=cumsum(E), y=cumsum(O))
abline(0,1, col="red")

# G test
G <- 2*nobs*sum(O*log(O/E),na.rm=T)
pval.G <- 1-pchisq(G,df=nmax-1)

# Chisq test == O teste est'a ruim pois 1 municipio possui 155 irregularidades em sa'ude.
chisq.test(x=O*nobs, p = E)

#removendo este municipio
chisq.test(x=O[1:nmax]*(nobs-1), p = E[1:nmax]/sum(E[1:nmax]))

# Calculando os erros esperados
pi00 = 903
pi01 = 50
pi11 = 252
pi10 = 46

p0 = pi00/(pi00+pi10)
p1 = pi11/(pi11+pi01)

r0 = 1-Eqn(p0,r,p)
r1 = 1-Eqn(p1,r,p)
r0
r1
1-r1-r0
