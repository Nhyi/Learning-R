par(mfrow=c(1,2))

N = 10000
y_rpois = rpois(N, lambda = 5)
hist(y_rpois,
     breaks = 100,
     xlab = "Number Rolled",
     main = "Poisson Distribution")

x = seq(0, 1, by = 0.01)
y = dbeta(x, 50, 50)
plot(x, y,
     main = "Beta Distribution",
     xlab = "Probability",
     ylab = "Frequency",)