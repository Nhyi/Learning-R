occ_simulation2 = function(n_year, lambda, alpha, beta, max_event){
  beta_events = matrix(rbeta(n_year * max_event, shape1 = alpha, shape2 = beta), nrow = n_year)
  n_events_per_year = rpois(n_year, lambda = lambda)
  for(i in which(n_events_per_year < max_event)) {
    beta_events[i, (n_events_per_year[i] + 1):max_event] = NA
  }
  cbind(1:n_year, beta_events)
}

occ_simulation2(n_year = 10, lambda = 10, alpha = 2, beta = 20, max_event = 5)