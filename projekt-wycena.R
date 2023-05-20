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


binomial_tree <- function(S_0, u, d_t, T){
  
  N <- T / d_t    # liczba kroków w drzewie, rozmiar macierzy
  A <- matrix(NA, N, N)
  for (i in 1:N){
    A[N:(N - i + 1), i] <- u**(seq(-(i - 1), (i - 1), 2))
  }
  S_T <- S_0 * A
  
  return(S_T)
}


S_T <- model(S_0, u, d, r, K, d_t, T)


wycena<-function(r,T,d,u,Vu,Vd,d_t){
  p <- (exp(r * d_t) - d) / (u - d)
  return(exp(-r * d_t) * (p *Vu + (1 - p) * Vd))
}


# opcja dla call@K

calle<-function(u,d,K,T,r,d_t,S_T){
  k <- T / d_t
  B <- matrix(NA, k, k)    # macierz payoff
  for (i in 1:k){
    B[i, k] <- max(S_T[i, k] - K, 0)    # ostatnia kolumna, czyli payoff dla S_T
  }

  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      B[j, i] <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1], d_t)
    }
  }
return(B)
  }

View(calle(u,d,K,T,r,d_t,S_T))



pute<-function(u,d,K,T,r,d_t,S_T){
  k <- T / d_t
  B <- matrix(NA, k, k)    # macierz payoff
  for (i in 1:k){
    B[i, k] <- max(K-S_T[i, k], 0)    # ostatnia kolumna, czyli payoff dla S_T
  }
  
  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      B[j, i] <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1], d_t)
    }
  }
  return(B)
}

european_option<-function(u,d,K,T,r,d_t,S_T, type='put'){
  k <- T / d_t
  B <- matrix(NA, k, k)    # macierz payoff
  for (i in 1:k){
    if(type=='put')
      B[i, k] <- max(K-S_T[i, k], 0)    # ostatnia kolumna, czyli payoff dla S_T
    else
      B[i, k] <- max(S_T[i, k]-K, 0)
  }
  
  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      B[j, i] <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1], d_t)
    }
  }
  return(B)
}

View(pute(u,d,K,T,r,d_t,S_T))


calla<-function(u,d,K,T,r,d_t,S_T){
  k <- T / d_t
  B <- matrix(NA, k, k)    # macierz payoff
  # for (i in 1:k){
  #  B[, k] <- pmax(S_T[, k]-K, 0) B[i, k] <- max(S_T[i, k] - K, 0)    # ostatnia kolumna, czyli payoff dla S_T
  # }
  # 
  B[, k] <- pmax(S_T[, k]-K, 0)
  
  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      a<-wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1], d_t)
      b<-max(S_T[j, i]-K,0)
      B[j, i] <- max(a,b)
    }
  }
  return(B)
}

View(calla(u,d,K,T,r,d_t,S_T))


puta<-function(u,d,K,T,r,d_t,S_T){
  k <- T / d_t
  B <- matrix(NA, k, k)    # macierz payoff
  for (i in 1:k){
    B[i, k] <- max(K-S_T[i, k], 0)    # ostatnia kolumna, czyli payoff dla S_T
  }
  
  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      a <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1])
      b<-max(K-S_T[j, i],0)
    B[j, i]<-max(a,b)
    }
  }
  return(B)
}

View(puta(u,d,K,T,r,d_t,S_T))


# badanie opłacalności? 

maksimum<-matrix(NA,k,k)

