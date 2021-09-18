#setting initial conditions
trials = 1000
lambda = 20

#creating empty list to store our variables
poisson_list = c()

#running the trials
for (i in 1:trials){
  #setting values for simulation of poisson
  X = 0
  P = 1
  
  #simulating poisson using runif
  while (P >= exp(-lambda)){
    U = runif(1, min = 0, max = 1)
    P = U * P
    X = X + 1
  }
  
  #appending values to the end of the list
  poisson_list = c(poisson_list, X)
}

#graphing our results
hist(poisson_list, main = "Histogram of Unif Pois", xlab = "Number")

#using inbuilt poisson function to compare against
trial = rpois(trials, lambda)

#plotting our results
hist(trial, main = "Histogram of RPois", xlab = "Number")