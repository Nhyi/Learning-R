library(ggplot2)

#setting initial conditions and number of trials
n = 1000000
one = two = three = 0


#iterating through a uniform distribution
for (i in 1:n){
  random = runif(1, min = 0, max = 1)

  #counting the number of times 1, 2 and 3 are rolled
  if (0.0 <= random && random < 0.4){
    one = one + 1
  } else if (0.4 <= random && random < 0.9){
    two = two + 1
  } else {
    three = three + 1
  }
  
}

#swapping to density
num_rolled = c((one/n), (two/n), (three/n))

#plotting results
par(mfrow = c(1,2))

#results using uniform distribution
barplot(num_rolled, main="Biased Dice Simulation (Uniform)", xlab="Number Rolled",  
        ylab="Density", names.arg=c("One","Two","Three"), 
        border="blue", col = heat.colors(3))

#comparison against inbuilt sample function
dice_probs = c(4/10, 5/10, 1/10)
trial = sample(1:3, size = n, replace = TRUE, prob = dice_probs)
expected = sample(1:3, size = n, replace = TRUE)
hist(trial, breaks = seq(0,3, 1), probability = TRUE, main = "Biased Dice Simulation (Sample)", 
     xlab = "Number Rolled", col = heat.colors(3))

