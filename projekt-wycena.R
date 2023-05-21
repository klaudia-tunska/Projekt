# trzymanie drzew w macierzy zaczynamy od dolnego lewego rogu, potem macierz payoff, funkcja pmax po wektorach
# macierz ceny, macierz payoff, macierz maksimów
# w europejskiej zerujemy payoffy do ostatniego, bo dopiero na końcu wykonujemy
# w amerykańskiej nie zerujemy payofffów
# robimy funkcje payoff dla call i put


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
  
  N <- T / d_t    # liczba kroków w drzewie, rozmiar macierzy
  A <- matrix(NA, N, N)
  for (i in 1:N){
    A[N:(N - i + 1), i] <- u**(seq(-(i - 1), (i - 1), 2))
  }
  S_T <- S_0 * A
  
  return(S_T)
}

S_T <- binomial_tree(S_0, u, d_t, T)

wycena <- function(u, d, r, d_t, V_u, V_d){
  
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p * V_u + (1 - p) * V_d))
}

european_option <- function(S_0, u, d, r, K, d_t, T, type = 'put'){
  
  N <- T / d_t    # liczba kroków w drzewie, rozmiar macierzy
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
  
  N <- T / d_t    # liczba kroków w drzewie, rozmiar macierzy
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

# Wykonanie opcji europejskiej put wynosi 6.2, natomiast amerykańskiej 6.37 (prawda).
# Wykonanie opcji europejskiej call wynosi 10, tak samo jak wykonanie takiej opcji amerykańskiej (prawda).

#Wizualizacja

pa<-american_option(S_0, u, d, r, K, d_t, T, type = 'put')
pheatmap(pa,Rowv = NA, border_color = "white")
heatmap(pa,Rowv = NA)
x<-pa[1:24,]
y<-seq(1,24,length=length(x))
plot(y,x)


# Zadanie 3
# Momenty wykonania opcji amerykańskich put
execution_time_put <-  which(american_option(S_0, u, d, r, K, d_t, T, type = 'put') == K - S_T, arr.ind = TRUE)
colnames(execution_time_put) <- c("moment", "czas")
execution_time_put <- execution_time_put[, 2:1]
execution_time_put

# Momenty wykonania opcji amerykańskich call
execution_time_call <-  which(american_option(S_0, u, d, r, K, d_t, T, type = 'call') == S_T - K, arr.ind = TRUE)
colnames(execution_time_call) <- c("moment", "czas")
execution_time_call <- execution_time_call[, 2:1]
execution_time_call


# Portfel

#portfel<-function(M){
#M<-american_option(S_0,u,d,r,K,d_t,T,type="call")
N <- T / d_t 
delta<-matrix(NA, N,N)
alfa<-matrix(NA, N,N)
  for (i in (N - 1):1){
    for (j in (N - i + 1):N){
      
  delta[j-1,i-1]<-(M[j,i]-M[j-1,i])/(S_T[j,i]-S_T[j-1,i]) #coś chyba tu nie do końca z wymiarem
  alfa[j-1,i-1]<-exp(-r*d_t)*(M[j,i]-delta[j-1,i-1]*S_T[j,i])
  
  }}
  
  
#  return(delta,alfa)
#}