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
  A <- matrix(NA, k, k)
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
