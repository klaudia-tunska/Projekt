library(pheatmap)

# nasze dane
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


binomial_tree <- function(S_0, u, d_t, T){
  
  N <- T / d_t +1  # liczba kroków w drzewie, rozmiar macierzy(dodajemy 1 gdyż R nie liczy od 0)
  A <- matrix(NA, N, N)
  for (i in 1:N){
    A[N:(N - i + 1), i] <- u**(seq(-(i - 1), (i - 1), 2))
  }
  S_T <- S_0 * A
  
  return(S_T)
}

S_T <- binomial_tree(S_0, u, d_t, T)
#View(S_T)

y<-S_T[,1:25]
x<-seq(1,25,length=length(y))
color<-c()
color[which(S_T>K)]<-"red"
color[which(S_T<K)]<-"blue"
plot(x,y, ylim=c(0,400),col=color, scale="log")# ta oś X żeby była od 0-2 jako czas w latach


wycena <- function(u, d, r, d_t, V_u, V_d){
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p * V_u + (1 - p) * V_d))
}

wycena<-Vectorize(wycena, c("V_u","V_d"))

european_option <- function(S_0, u, d, r, K, d_t, T, type = 'put'){
  
  N <- T / d_t+1    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)    # macierz w momencie S_t
  B <- matrix(0, N, N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == 'put')
    B[, N] <- pmax(K - S_T[, N], 0)
  else
    B[, N] <- pmax(S_T[, N] - K, 0)
  
  for (i in (N - 1):1){
    for (j in (N - i + 1):N){
      B[j, i] <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
    }
  }
  
  return(B)
}


EU_put <- round(european_option(S_0, u, d, r, K, d_t, T, type = 'put')[T/d_t+1, 1], 2)    # zadanie 1
EU_call <- round(european_option(S_0, u, d, r, K, d_t, T, type = 'call')[T/d_t+1, 1], 2)

#View(european_option(S_0, u, d, r, K, d_t, T, type = 'put'))

american_option <- function(S_0, u, d, r, K, d_t, T, type = 'put'){
  
  N <- T / d_t+1   # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)    # macierz w momencie S_t
  B <- matrix(0, N, N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == 'put'){
    B[, N] <- pmax(K - S_T[, N], 0)
    
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(K - S_T[j, i], 0)
        B[j, i] <- max(a, b)
      }
    }
  }
  else{
    B[, N] <- pmax(S_T[, N] - K, 0)
  
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(S_T[j, i] - K, 0)
        B[j, i] <- max(a, b)
      }
    }
  }
    return(B)
  }


AM_put <- round(american_option(S_0, u, d, r, K, d_t, T, type = 'put')[T/d_t+1, 1], 2)    # zadanie 2
AM_call <- round(american_option(S_0, u, d, r, K, d_t, T, type = 'call')[T/d_t+1, 1], 2)

# Cena opcji europejskiej put wynosi 6.2, natomiast amerykańskiej 6.37 (prawda).
# Cena opcji europejskiej call wynosi 10, tak samo jak wykonanie takiej opcji amerykańskiej (prawda).

#####
ec<-european_option(S_0, u, d, r, K, d_t, T, type = 'call')
ep<-european_option(S_0, u, d, r, K, d_t, T, type = 'put')
ac<-american_option(S_0, u, d, r, K, d_t, T, type = 'call')
ap<-american_option(S_0, u, d, r, K, d_t, T, type = 'put')
View(ec)
View(ep)
View(ac)
View(ap)
#####
#Wizualizacja

pe<-european_option(S_0, u, d, r, K, d_t, T, type = 'call')
pa<-american_option(S_0, u, d, r, K, d_t, T, type = 'call')
pheatmap(pa,Rowv = NA, Colv=NA,border_color = "white",show_rownames=TRUE,show_colnames=TRUE)
library(RColorBrewer)

heatmap(pa,Rowv = NA,Colv = NA, scale="column", xlab="something",
        ylab="", main="A title",col= colorRampPalette(brewer.pal(8, "Oranges"))(25))

y<-pa[1:25,]
x<-seq(1,25,length=length(y))
plot(x,y, ylim=c(0,59), col=c(9,0,1), type="o")


# Zadanie 3
# Momenty wykonania opcji amerykańskich put
S_T <- binomial_tree(S_0, u, d_t, T)
execution_time_put <-  which(american_option(S_0, u, d, r, K, d_t, T, type = 'put') == K - S_T, arr.ind = TRUE)
colnames(execution_time_put) <- c("moment", "czas")
execution_time_put <- execution_time_put[, 2:1]
execution_time_put

# Momenty wykonania opcji amerykańskich call
execution_time_call <-  which(american_option(S_0, u, d, r, K, d_t, T, type = 'call') == S_T - K, arr.ind = TRUE)
colnames(execution_time_call) <- c("moment", "czas")
execution_time_call <- execution_time_call[, 2:1]
execution_time_call


moments<-american_option(S_0, u, d, r, K, d_t, T, type = 'put') == K - S_T
moments[moments==TRUE]<-1
moments[moments==FALSE]<-0
heatmap(moments,Rowv = NA,Colv = NA) #fajnie byłoby jednym kolorem pozaznaczać te miejsca i jakoś obrócić ta mape??
View(moments)

# Zadanie 4

# Wrażliwość na zmianę ceny wykonania K
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
T <- 2
#gdzieś wcześniej są modyfikacje do funkcji, więc wczyrujemy jeszcze raz "czystą" funkcję 
european_option <- function(S_0, u, d, r, K, d_t, T, type){
  
  N <- T / d_t+1    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)    # macierz w momencie S_t
  B <- matrix(0, N, N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == 'put')
    B[, N] <- pmax(K - S_T[, N], 0)
  else
    B[, N] <- pmax(S_T[, N] - K, 0)
  
  for (i in (N - 1):1){
    for (j in (N - i + 1):N){
      B[j, i] <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
    }
  }
  
  return(B)
} 
american_option <- function(S_0, u, d, r, K, d_t, T, type){
  
  N <- T / d_t+1   # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)    # macierz w momencie S_t
  B <- matrix(0, N, N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == 'put'){
    B[, N] <- pmax(K - S_T[, N], 0)
    
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(K - S_T[j, i], 0)
        B[j, i] <- max(a, b)
      }
    }
  }
  else{
    B[, N] <- pmax(S_T[, N] - K, 0)
    
    for (i in (N - 1):1){
      for (j in (N - i + 1):N){
        a <- wycena(u, d, r, d_t, B[j - 1, i + 1], B[j, i + 1])
        b <- max(S_T[j, i] - K, 0)
        B[j, i] <- max(a, b)
      }
    }
  }
  return(B)
}

european_option_K<-Vectorize(european_option,"K")
american_option_K<-Vectorize(american_option, "K")

K<-40:50
ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()

for (i in 1:length(K)){
  ceny_pe[i]<-european_option_K(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_ce[i]<-european_option_K(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_pa[i]<-american_option_K(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ca[i]<-american_option_K(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
}



ceny_pe
ceny_ce
ceny_pa
ceny_ca


plot(K,ceny_pe,col="green", ylab = "Wartość", xlab = "K", main="Zależność ceny opcji od , K")
lines(K,ceny_ce, col="blue", type="p")
lines(K,ceny_pa,col="black")
lines(K,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))





# Wrażliwość na zmianę zapadalności T
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


european_option_T<-Vectorize(european_option,"T")
american_option_T<-Vectorize(american_option, "T")


T<-seq(1,100,by=1)
ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()


for (i in 1:length(T)){
  ceny_pe[i]<-european_option_T(S_0, u, d, r, K, d_t, T, type = 'put')[[i]][T[i]/d_t+1,1]
  ceny_ce[i]<-european_option_T(S_0, u, d, r, K, d_t, T, type = 'call')[[i]][T[i]/d_t+1,1]
  ceny_pa[i]<-american_option_T(S_0, u, d, r, K, d_t, T, type = 'put')[[i]][T[i]/d_t+1,1]
  ceny_ca[i]<-american_option_T(S_0, u, d, r, K, d_t, T, type = 'call')[[i]][T[i]/d_t+1,1]
}



ceny_pe
ceny_ce
ceny_pa
ceny_ca


plot(T,ceny_pe,col="green",ylim=c(0,20), ylab = "Wartość", xlab = "Zapadalnoś T", main="Zależność ceny opcji od zapadalności, T")
lines(T,ceny_ce, col="blue", type="p")
lines(T,ceny_pa,col="black")
lines(T,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))


plot(maturity,eu_call_T,type="l",col=2,xlab="Zapadalność [w latach]",ylab="Cena opcji",
     main="Zależność ceny opcji od zapadalności")
lines(maturity,eu_put_T,type="l",col=3)
lines(maturity,am_call_T,type="l",col=4)
lines(maturity,am_put_T,type="l",col=5)
legend(0,2.75,legend=c("eu call","eu put","am call","am put"),lty=3,pch=18,col=c(2,3,4,5),cex=0.7)
?legend
#legenda nie działa?????



# Wrażliwość na zmianę S_0
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


european_option_S0<-Vectorize(european_option,"S_0")
american_option_S0<-Vectorize(american_option, "S_0")


S_0<-40:60
ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()


for (i in 1:length(S_0)){
  ceny_pe[i]<-european_option_S0(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ce[i]<-european_option_S0(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_pa[i]<-american_option_S0(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ca[i]<-american_option_S0(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
}


ceny_pe
ceny_ce
ceny_pa
ceny_ca

plot(S_0,ceny_pe,col="green",ylim=c(0,20), ylab = "Wartość",
     xlab = "Wartość S0", main="Zależność ceny opcji od S0")
lines(S_0,ceny_ce, col="blue", type="p")
lines(S_0,ceny_pa,col="black")
lines(S_0,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))



# Wrażliwość na zmianę sigma
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


european_option_sigma<-Vectorize(european_option,c("u","d"))
american_option_sigma<-Vectorize(american_option,c("u","d"))

sigma<-seq(0.05,0.5,by=0.05)
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))

ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()


for (i in 1:length(sigma)){
  ceny_pe[i]<-european_option_sigma(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ce[i]<-european_option_sigma(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_pa[i]<-american_option_sigma(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ca[i]<-american_option_sigma(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
}


ceny_pe
ceny_ce
ceny_pa
ceny_ca


plot(sigma,ceny_pe,col="green",ylim=c(0,20), ylab = "Wartość",
     xlab = "Wartość sigmy", main="Zależność od sigmy")
lines(sigma,ceny_ce, col="blue", type="p")
lines(sigma,ceny_pa,col="black")
lines(sigma,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))



# Wrażliwość na zmianę r
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


european_option_r<-Vectorize(european_option,"r")
american_option_r<-Vectorize(american_option,"r")

r<-seq(0.01,0.2,by=0.01)

ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()


for (i in 1:length(r)){
  ceny_pe[i]<-european_option_r(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ce[i]<-european_option_r(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_pa[i]<-american_option_r(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ca[i]<-american_option_r(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
}


ceny_pe
ceny_ce
ceny_pa
ceny_ca


plot(r,ceny_pe,col="green",ylim=c(0,20), ylab = "Wartość", xlab = "Wartość r", main="Zależność od r")
lines(r,ceny_ce, col="blue", type="p")
lines(r,ceny_pa,col="black")
lines(r,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                        "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))


# Zadanie 5

# Wrażliwość na d_t
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


european_option_dt<-Vectorize(european_option,  c("u","d","d_t"))
american_option_dt<-Vectorize(american_option, c("u","d","d_t"))

T<-1
d_t<-seq(0.01,1, by=0.04)
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))

ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()


for (i in 1:length(d_t)){
ceny_pe[i]<-european_option_dt(S_0, u, d, r, K, d_t, T, type = 'put')[[i]][T/d_t[i]+1,1]
ceny_ce[i]<-european_option_dt(S_0, u, d, r, K, d_t, T, type = 'call')[[i]][T/d_t[i]+1,1]
ceny_pa[i]<-american_option_dt(S_0, u, d, r, K, d_t, T, type = 'put')[[i]][T/d_t[i]+1,1]
ceny_ca[i]<-american_option_dt(S_0, u, d, r, K, d_t, T, type = 'call')[[i]][T/d_t[i]+1,1]
}

ceny_pe
ceny_ce
ceny_pa
ceny_ca

plot(d_t,ceny_pe,col="green",ylim=c(0,10), ylab = "Wartość", xlab = "Krok", main="Zależność ceny opcji od d_t")
lines(d_t,ceny_ce, col="blue", type="p")
lines(d_t,ceny_pa,col="red")
lines(d_t,ceny_ca, col="magenta")
legend("bottomright", c("Europejska put","Europejska call","Amerykańska put",
                        "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","red","magenta"))


# Zadanie 6

d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2

portfel<-function(M, S_T){
N <- T / d_t + 1
delta<-matrix(NA, N,N)
alfa<-matrix(NA, N,N)
  for (i in 1:N-1){ 
    for (j in N:2){
  delta[j,i]<-(M[j-1,i+1]-M[j,i+1])/(S_T[j-1,i+1]-S_T[j,i+1])
  alfa[j,i]<-exp(-r*d_t)*(M[j-1,i+1]-delta[j,i]*S_T[j-1,i+1])
  
  }}
  x<-list(delta,alfa)
  return(x)
}

S_T <- binomial_tree(S_0, u, d_t, T)

pe<-european_option(S_0,u,d,r,K,d_t,T,type="put")
ce<-european_option(S_0,u,d,r,K,d_t,T,type="call")
pa<-american_option(S_0,u,d,r,K,d_t,T,type="put")
ca<-american_option(S_0,u,d,r,K,d_t,T,type="call")

portfel(pe, S_T)
portfel(ce, S_T)
portfel(pa, S_T)
portfel(ca, S_T)
