dice_function <- function(probs, xvals, trials){
  
  data = NULL
  cumulative = cumsum(probs)
  
  for(i in 1:trials){
    
    uniform = runif(1)
    data = c(data, findInterval(uniform, cumulative))
    ylim = max(probs)*trials
    
  }
  
  pdf (file = "D:\\Documents\\Programming\\R Projects\\Shree's Project\\Week 4\\Dice Plot (100000 Trials).pdf", width = 4, height = 4)
  
  barplot(table(data),
          col = rgb(0.8, 0.1, 0.1, 0.6),
          xlab = "X Values",
          ylab = "Frequency",
          main = "Histogram of Manual Dice Rolling",
          ylim = c(0, ylim))
  
  theoretical = sample((min(xvals)):max(xvals), size = trials, replace = TRUE, prob = probs)
  
  hist(theoretical, col = rgb(0.8, 0.1, 0.1, 0.6),
       main = "Histogram of Theoretical Dice Rolling", xlab = "X Values")
  
  dev.off()

}

poisson_function <- function(lambda, trials){
  
  data = NULL
  
  for (i in 1:trials){
    
    X = 0
    P = 1
    
    while (P >= exp(-lambda)){
      
      U = runif(1, min = 0, max = 1)
      P = U * P
      X = X + 1
      
    }
    
    data = c(data, X)
    
  }
  
  pdf (file = "D:\\Documents\\Programming\\R Projects\\Shree's Project\\Week 4\\Poisson Plot (100000 Trials).pdf", width = 4, height = 4)
  
  barplot(table(data),
          col = rgb(0.8, 0.1, 0.1, 0.6),
          xlab = "X Values",
          ylab = "Frequency",
          main = "Histogram of Manual Poisson")
  
  theoretical = rpois(trials, lambda)
  hist(theoretical, main = "Histogram of Theoretical Poisson",
       col = rgb(0.8, 0.1, 0.1, 0.6), xlab = "X Values", ylab = "Frequency", breaks = 30)
  
  dev.off()
}

a = c(0.05, 0.05, 0.4, 0.3, 0.2)
b = c(0, 1, 2, 3, 4)
manual_dice = dice_function(a, b, 100000)

manual_poisson = poisson_function(20, 100000)