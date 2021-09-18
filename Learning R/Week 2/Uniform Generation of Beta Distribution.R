#setting initial conditions
N = 10000
alpha = 5
beta = 5

par(mfrow = c(1,2))

#running the simulation
U = runif(N)
b_rand = qbeta(U, alpha, beta)

#plotting the results
hist(b_rand, xlab = "Value", col="skyblue", main = "Runif For Beta Distribution")

#draw N beta distributed values
y_rbeta = rbeta(N, shape1 = alpha, shape2 = beta)

#plot of randomly drawn beta density
plot(density(y_rbeta), main = "Beta Distribution With Rbeta")