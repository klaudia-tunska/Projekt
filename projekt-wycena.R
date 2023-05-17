# trzymanie drzew w macierzy zaczynamy od dolnego lewego rogu, potem macierz payoff, funkcja pmax po wektorach
# macierz ceny, macierz payoff, macierz maksimów
# w europejskiej zerujemy payoffy do ostatniego, bo dopiero na końcu wykonujemy
# w amerykańskiej nie zerujemy payofffów
# robimy funkcje payoff dla call i put


# nasze dane
d_t <- 1 / 12
sigma <- 0.3
u <- exp(sigma * sqrt(d_t))
d <- exp(-sigma * sqrt(d_t))
S_0 <- 50
r <- 0.02
K <- 48
T <- 2


model <- function(S_0, u, d, r, K, d_t, T){
  
  k <- T / d_t    # liczba kroków w drzewie, rozmiar macierzy
  A <- matrix(0, k, k)
  A[k, 1] <- 1
  for (i in 2:k){
    for (j in (k - i + 1):k){
      if (j == k) {
        A[j, i] <- d * A[j, i - 1]    # krok w dół (ale stosujemy tylko przy ostatniej wartości w momencie t)
      }
      else {
        A[j, i] <- u * A[j + 1, i - 1]    # krok w górę
      }
    }
  }
  S_T = S_0 * A
  
  return(S_T)


}


S_T <- model(S_0, u, d, r, K, d_t, T)


wycena<-function(r,T,d,u,Vu,Vd){
  p <- (exp(r * T) - d) / (u - d)
  return(exp(-r * d_t) * (p *Vu + (1 - p) * Vd))
}


# opcja dla call@K

calle<-function(u,d,K,T,r,d_t,S_T){
  k <- T / d_t
  B <- matrix(0, k, k)    # macierz payoff
  for (i in 1:k){
    B[i, k] <- max(S_T[i, k] - K, 0)    # ostatnia kolumna, czyli payoff dla S_T
  }

  for (i in (k - 1):1){
    for (j in (k - i + 1):k){
      B[j, i] <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1])
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
      B[j, i] <- wycena(r,T,d,u,B[j - 1, i + 1], B[j, i + 1])
    }
    }
  return(B)
}

View(pute(u,d,K,T,r,d_t,S_T))
