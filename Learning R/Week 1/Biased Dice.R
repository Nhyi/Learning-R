dice_probs = c(4/10, 5/10, 1/10)
trial = sample(1:3, size = 1000000, replace = TRUE, prob = dice_probs)
par(mfrow = c(1,2))
hist(trial, breaks = seq(0,3, 1), probability = TRUE, col = rainbow(6),
     main = "Results of Biased Dice Rolls", xlab = "Number")

expected = sample(1:3, size = 1000000, replace = TRUE)
hist(expected, breaks = seq(0,3, 1), probability = TRUE, col = rainbow(6),
     main = "Results of Unbiased Dice Rolls", xlab = "Number")