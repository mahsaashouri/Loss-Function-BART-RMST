
## Function calculating the probability of each tree structure; 
## d: row number in Dmat matrix (matrix involving all the 26 possible structure)

DProb <- function(alpha, beta, d){
  if(d == 1)
    prob <- 1-alpha
  else{
    if(d == 2){
      m = r = h = f = 0; l = 2;
    }else if(d == 3 | d == 4){
      m = l = 1; r = f = 0; h = 2;
    }else if(d == 5){
      m = 2; l = r = f = 0; h = 4;
    }else if(d == 6 | d == 7 | d == 9 | d == 10){
      m = l = r = h = 1; f = 2;
    }else if(d == 8 | d == 11){
      m = l = 1; r = 2; h = 0; f = 4;
    }else if(d == 12 | d == 13 | d == 14 | d == 15){
      m = 2; l = 0; r = 1; h = 3; f = 2;
    }else if(d == 16 | d == 17 | d == 18 | d == 19 | d == 20 | d == 21){
      m = r = h = 2; l = 0; f = 4;
    }else if(d == 22 | d == 23 | d == 24 | d == 25){
      m = 2; l = 0; r = 3; h = 1; f = 6;
    }else{
      m = 2; l = 0; r = 4; h = 0; f = 8;
    }
    prob <- alpha*(alpha/(2^beta))^m*(1-(alpha/(2^beta)))^l*(alpha/(3^beta))^r*(1-(alpha/(3^beta)))^h*(1-(alpha/(4^beta)))^f
  }   
  return(prob)
}


## test the sum of probabilities
SumprobD <- c()
for(i in 1:26){
  SumprobD[i] <- DProb(alpha = .95, beta = 2, i)
  
}
sum(SumprobD)
