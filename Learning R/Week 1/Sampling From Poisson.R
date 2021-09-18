for (i in 1:10){
  n = rpois (1,lambda = 2);
  print(n);
  
    if (n > 0) {
      for (i in 1:n){
        b = rbeta (1, shape1 = 1, shape2 = 1.5);
        print(b)
      }
    }
}