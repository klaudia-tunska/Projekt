# trzymanie drzew w macierzy zaczynamy od dolnego lewego rogu, potem macierz payoff, funkcja pmax po wektorach
# macierz ceny, macierz payoff, macierz maksimów
# w europejskiej zerujemy payoffy do ostatniego, bo dopiero na końcu wykonujemy
# w amerykańskiej nie zerujemy payofffów
# robimy funkcje payoff dla call i put


dt<-1/12
sigma<-0.3
u<-exp(sigma*sqrt(dt))
d<-exp(-sigma*sqrt(dt))
so<-50
r<-0.02
k<-48
t<-2
