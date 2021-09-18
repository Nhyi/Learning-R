library(tidyverse)

occ_simulation <- function(nyears, lambda, alpha, beta){
  
  data_matrix <- matrix(, nrow = nyears, ncol = 6)
  
  for (z in 1:nyears){
    data_matrix[z][1] <- z
  }
  
  for (yr in 1:nyears){
    
    poisson_sim = rpois(1, lambda)
    
    for (number_of_events in poisson_sim){
      
      if (number_of_events == 1){
        
        beta_sim = rbeta(1, alpha, beta)
        data_matrix[yr, 2] <- beta_sim
        
      } else if (number_of_events == 2){
        
          for (i in 2:3){
            
            beta_sim = rbeta(1, alpha, beta)
            data_matrix[yr, i] <- beta_sim
            
          }
        
      } else if (number_of_events == 3){
          
          for (i in 2:4){
            
            beta_sim = rbeta(1, alpha, beta)
            data_matrix[yr, i] <- beta_sim
        
          }
        
      } else if (number_of_events == 4){
        
          for (i in 2:5){
            
            beta_sim = rbeta(1, alpha, beta)
            data_matrix[yr, i] <- beta_sim
        
          }
          
      } else{
        
          for (i in 2:6){
            
            beta_sim = rbeta(1, alpha, beta)
            data_matrix[yr, i] <- beta_sim
        }
    
      }
      
    }
    
  }
  
  sorted_matrix <- cbind(data_matrix[,1],t(apply(data_matrix[,2:6],1,function(x) sort(x))))
  
G <- sorted_matrix %>% as.data.frame %>%
  pivot_longer(-V1) %>%
  ggplot(aes(x=factor(V1),y=value,color=name,group=name))+
  geom_point()+
  labs(color='Column',x='Time (Years)', y ='Probability')+
  theme_bw()
return(G)
}

manual = occ_simulation(10, 10, 2, 20)
manual