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


a<-matrix(NA, 24, 24) 
a[24,1] #lewy dolny róg

for (n in 23:1){
a[n,25-n]<-u^(24-n)
}
a
k<-23:1
k
a[24,25-k]<-d^(24-k)
a[24,]
a


a[w,c+1]<-a[w+1,c]*a[w,c]
a[24,2]*a[23,2]
c=3
w=23
a[24,c-1]*a[w,c-1]

n<-1
for (c in 2:23){
  for (w in 1:n){
    a[24-w,c+1]<-a[24+1-w,c]*a[24-w,c]
  }
  n<-n+1
} 

View(a)

payoff_call<-function(st,k){
  max(st-k,0)
}

pec<-matrix(NA,24,24)
for (i in 1:24){
pec[i,i]<-max(a[i,i]-k,0)
}


k <- 10 # rozmiar macierzy (później = 24)
A <- matrix(0, k, k)
A[k, 1] <- 1
for (i in 2:k){
  for (j in (k - i + 1):k){
    if (j == k) {
      A[j, i] <- d * A[j, i - 1]
    }
    else {
      A[j, i] <- u * A[j + 1, i - 1]
    }
  }
}
View(A)
