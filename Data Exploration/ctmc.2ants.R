######################################
##
##  CTMC Coding
##
######################################


#Goal: simulate interaction chain between two ants

Ti = 10 #number of interactions
ant =3 #number of states

n = 1 + #no interactions
  choose(ant, 2) + #two way interactions
  
  choose(ant, 3)  #three way interactions

#homes for Embedded Chain and Rate Transitions

EC = rep(NA, Ti)
RT = rep(NA, Ti)

#transitions rates (not probabilities!)
alpha = matrix(1, n, n, Ti) #(1->2, 2->1)

for(i in 1:n){
  alpha[i, i] = 0
}

v= rep(NA, n) #rate for state i

for(i in 1:n){
  v[i] = sum(alpha[i, ])
}

P = matrix(NA, n, n) #transition probability matrix

for(i in 1:n){
  for(k in 1:n){
    P[i, k] = alpha[i, k] / v[i]
  }
}


#initialize EC
EC[1] = sample(1:n, 1, prob = rep(1/n, n))
RT[1] = rexp(1, rate = v[EC[1]])

for(t in 2:Ti){
  EC[t] = sample(1:n, 1, prob = P[EC[t-1],])
  RT[t] = rexp(1, rate = v[EC[t]])
}


##plot of interactions

  
  t = cumsum(RT)  
  X = EC

  plot(0,0,xlab="t",ylab="State",ylim=c(1,max(X)),xlim=c(0,max(t)),type="n")
      lines(x = c(0, t[1]), y = c(X[1], X[1]))
  
   for(i in 2:Ti){
    lines(x = c(t[i-1], t[i]), y = c(X[i], X[i]))
   }

  
##estimate alpha values
  
  m.ij = matrix(NA, n, n) #matrix of number of times chain moves from state i to state j

  
