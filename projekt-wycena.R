# Konstruujemy drzewo niezależnie od ścieżki (PATH-INDEPENDENT) jako macierz, 
# gdzie t-ta kolumna to S_t moment wykonania akcji

# Dane
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2

matrix_size <- function(d_t, T){    # rozmiar macierzy kwadratowej N, gdzie uwzględniamy moment S_0
  return(T / d_t + 1)
}

N <- matrix_size(d_t, T)

binomial_tree <- function(S_0, u, d_t, T){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  A <- matrix(NA, nrow = N, ncol = N)
  
  for (i in 1:N){
    A[N:(N - i + 1), i] <- u**(seq(-(i - 1), (i - 1), 2))
  }
  
  return(S_0 * A)
}

S_T <- binomial_tree(S_0, u, d_t, T)
#View(S_T)

# Path-independent tree
y <- S_T[, 1:25]
x <- seq(1, 25, length = length(y))
colors <- ifelse(S_T > K, "#DEA435", "#40DE8E")
plot(x, y, log = "y", type = "p", col = colors, pch = 16, xlab = "Moments", ylab = "", main = "Path-independent tree", yaxt = "n")
legend("topleft", legend = c("S_T > K", "S_T < K"), col = c("#DEA435", "#40DE8E"), pch = 16, cex = 1.1, bty = "n")

wycena <- function(u, d, r, d_t, V_u, V_d){
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p * V_u + (1 - p) * V_d))
  }
wycena<-Vectorize(wycena, c("V_u","V_d"))

european_option <- function(S_0, u, d, r, K, d_t, T, type = "put"){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)
  B <- matrix(0, nrow = N, ncol = N)    # macierz payoff
  B[is.na(S_T)] <- NA
  ifelse(type == "put",
         B[, N] <- pmax(K - S_T[, N], 0),    # opcja put
         B[, N] <- pmax(S_T[, N] - K, 0))    # opcja call

  for (i in (N - 1):1){
    B[(N - i + 1):N, i] <- wycena(u, d, r, d_t, B[(N - i):(N - 1), i + 1], B[(N - i + 1):N, i + 1])
  }

  return(B)
}

#View(european_option(S_0, u, d, r, K, d_t, T, type = 'put'))

american_option <- function(S_0, u, d, r, K, d_t, T, type = "put"){
  
  N <- matrix_size(d_t, T)    # liczba kroków w drzewie, rozmiar macierzy
  S_T <- binomial_tree(S_0, u, d_t, T)
  B <- matrix(0, nrow = N, ncol = N)    # macierz payoff
  B[is.na(S_T)] <- NA
  if(type == "put"){
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

K<-20:80
ceny_pe<-c()
ceny_ce<-c()
ceny_pa<-c()
ceny_ca<-c()

for (i in 1:length(K)){
  ceny_pe[i]<-european_option_K(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ce[i]<-european_option_K(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
  ceny_pa[i]<-american_option_K(S_0, u, d, r, K, d_t, T, type = 'put')[25,i]
  ceny_ca[i]<-american_option_K(S_0, u, d, r, K, d_t, T, type = 'call')[25,i]
}



ceny_pe
ceny_ce
ceny_pa
ceny_ca


plot(K,ceny_pe,col="green", ylim=c(0,50) ,ylab = "Cena", xlab = "Cena spot",
     main="Zależność ceny opcji od ceny wykonania, K")
lines(K,ceny_ce, col="blue", pch=4, type = "p")
lines(K,ceny_pa,col="black")
lines(K,ceny_ca, col="red")
legend("topright", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","x","-","-"), col=c("green","blue","black","red"))





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


T<-seq(0,40,by=1)
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


plot(T,ceny_pe,col="green",ylim=c(0,20), ylab = "Cena", xlab = "Zapadalnoś T", main="Zależność ceny opcji od zapadalności, T")
lines(T,ceny_ce, col="blue", type="p")
lines(T,ceny_pa,col="black")
lines(T,ceny_ca, col="red")
legend("topleft", c("Europejska put","Europejska call","Amerykańska put",
                    "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","black","red"))



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

plot(S_0,ceny_pe,col="green",ylim=c(0,20), ylab = "Cena",
     xlab = "Wartość ceny spot", main="Zależność ceny opcji od ceny spot")
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

sigma<-seq(0.1,5,by=0.5)
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


plot(sigma,ceny_pe,col="green",ylim=c(0,50), ylab = "Cena",
     xlab = "Wartość sigmy", main="Zależność od sigmy")
lines(sigma,ceny_ce, col="blue", type="p")
lines(sigma,ceny_pa,col="black")
lines(sigma,ceny_ca, col="red")
legend("bottomright", c("Europejska put","Europejska call","Amerykańska put",
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

r<-seq(0,1,by=0.05)

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


plot(r,ceny_pe,col="green",ylim=c(0,50), ylab = "Cena", xlab = "Wartość r", main="Zależność ceny od r")
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

plot(d_t,ceny_pe,col="green",ylim=c(2,9), ylab = "Cena", xlab = "Krok", main="Zależność ceny opcji od liczby kroków")
lines(d_t,ceny_ce, col="blue", type="p")
lines(d_t,ceny_pa,col="red")
lines(d_t,ceny_ca, col="magenta")
legend("bottomright", c("Europejska put","Europejska call","Amerykańska put",
                        "Amerykańska call"),pch=c("o","o","-","-"), col=c("green","blue","red","magenta"))


# Zadanie 6

d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))
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

portfel_pe<-portfel(pe, S_T)
portfel_ce<-portfel(ce, S_T)
portfel_pa<-portfel(pa, S_T)
portfel_ca<-portfel(ca, S_T)

delta_pe<-portfel_pe[[1]][1:25,]
delta_ce<-portfel_ce[[1]][1:25,]
delta_pa<-portfel_pa[[1]][1:25,]
delta_ca<-portfel_ca[[1]][1:25,]

alfa_pe<-portfel_pe[[2]][1:25,]
alfa_ce<-portfel_ce[[2]][1:25,]
alfa_pa<-portfel_pa[[2]][1:25,]
alfa_ca<-portfel_ca[[2]][1:25,]

x<-seq(1,25,length=625)

par(mfrow = c(1, 2))    
plot(x,delta_pe, col=11,pch=19, main = "Wartości delty dla europejskiej opcji put",
     ylab = "Wartość delty", xlab = "Krok")
plot(x,alfa_pe,col=12,pch=19, main = "Wartości alfy dla europejskiej opcji put",
     ylab = "Wartość alfy", xlab = "Krok")

par(mfrow = c(1, 2))    
plot(x,delta_ce,col=11,pch=19, main = "Wartości delty dla europejskiej opcji call",
     ylab = "Wartość delty", xlab = "Krok")
plot(x,alfa_ce,col=12,pch=19, main = "Wartości alfy dla europejskiej opcji call",
     ylab = "Wartość alfy", xlab = "Krok")

par(mfrow = c(1, 2))    
plot(x,delta_pa,col=11,pch=19, main = "Wartości delty dla amerykańskiej opcji put",
     ylab = "Wartość delty",  xlab = "Krok")
plot(x,alfa_pa,col=12,pch=19, main = "Wartości alfy dla amerykańskiej opcji put",
     ylab = "Wartość alfy", xlab = "Krok")

par(mfrow = c(1, 2))    
plot(x,delta_ca,col=11,pch=19, main = "Wartości delty dla amerykańskiej opcji call",
     ylab = "Wartość delty", xlab = "Krok")
plot(x,alfa_ca,col=12,pch=19, main = "Wartości alfy dla amerykańskiej opcji call",
     ylab = "Wartość alfy", xlab = "Krok")
