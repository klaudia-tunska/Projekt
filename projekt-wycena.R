# nasze dane
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))    # d = 1 / u = u**(-1)
S_0 <- 50
r <- 0.02
K <- 48
T <- 2

library(pheatmap)

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
View(S_T)

wycena <- function(u, d, r, d_t, V_u, V_d){
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p * V_u + (1 - p) * V_d))
}


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

EU_put <- round(european_option(S_0, u, d, r, K, d_t, T, type = 'put')[T / d_t, 1], 2)    # zadanie 1
EU_call <- round(european_option(S_0, u, d, r, K, d_t, T, type = 'call')[T / d_t, 1], 2)


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


AM_put <- round(american_option(S_0, u, d, r, K, d_t, T, type = 'put')[T / d_t, 1], 2)    # zadanie 2
AM_call <- round(american_option(S_0, u, d, r, K, d_t, T, type = 'call')[T / d_t, 1], 2)

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

y<-pe[1:24,]
x<-seq(1,24,length=length(y))
plot(x,y, ylim=c(0,59))
lines(1:24,rep(K,time=24), type = "l")


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


# Wrażliwość na zmianę zapadalności T



# Zadanie 5

# Wrażliwość na d_t
interval <- seq(12,365,1)



# Zadanie 6

portfel<-function(M){
N <- T / d_t 
delta<-matrix(NA, N,N)
alfa<-matrix(NA, N,N)
  for (i in (N):1){ # czy nie powinno być do 2 bo gdy i=1 to wpisujemy potem w zerową kolumnę?
    for (j in (N - i + 1):N){
      
  delta[j-1,i-1]<-(M[j,i]-M[j-1,i])/(S_T[j,i]-S_T[j-1,i])
  alfa[j-1,i-1]<-exp(-r*d_t)*(M[j,i]-delta[j-1,i-1]*S_T[j,i])
  
  }}
  
  x<-list(delta,alfa)
  return(x)
}

M<-american_option(S_0,u,d,r,K,d_t,T,type="call")
portfel(M)
#może by zobrazoawać jakoś tą liste na wykresie?
########inaczej numeracja w portfelu, bo czy tamta ok była?
portfel<-function(M){
  N <- T / d_t 
  delta<-matrix(NA, N,N)
  alfa<-matrix(NA, N,N)
  for (i in 1:(N-1)){
    for (j in (N - i + 1):N){
      
      delta[j,i]<-(M[j-1,i+1]-M[j,i+1])/(S_T[j-1,i+1]-S_T[j,i+1])
      alfa[j,i]<-exp(-r*d_t)*(M[j-1,i+1]-delta[j,i]*S_T[j-1,i+1])
      
    }}
  
  x<-list(delta,alfa)
  return(x)
}

M<-american_option(S_0,u,d,r,K,d_t,T,type="call")
portfel(M)

