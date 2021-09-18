# Write a function that samples from an n-dimensional vector of xvals(n) and probs(n), m times, only using the uniform distribution. 
# 
# This function applies to any discrete probability vector that I feed into (which can be derived from a coin, a poisson, 
#or even a continuous pdf like beta (which requires numerical integration)).


sampling_function <- function(m, xvals, probs){

  uniform = runif(m, min = min(xvals), max = (xvals))
  
}

xvalues = c(1, 2, 3, 4, 5)
probs = c(0.2, 0.2, 0.2, 0.2, 0.2)

sample_result = do.call(sampling_function, xvalues, probs)