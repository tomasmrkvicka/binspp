dgenp=function (x, lambda, theta){
  norm=1
  if(lambda<0){
    n=ceiling(-theta/lambda)
    norm=0
    for(k in 0:n) {
      norm =  norm+
        max(0,theta*(theta + lambda*k)^(k-1)*exp(-theta-lambda*k) /factorial(k)) }
  }
  if((lambda > 1) || (lambda < -1) || (theta<0)){
    stop("Parameters not allowed")}
  return(max(0,theta*(theta + lambda*x)^(x-1)*exp(-theta-lambda*x) /factorial(x))/norm)
}


rgenp=function (n, lambda, theta)
{
  if((lambda > 1) || (lambda < -1) || (theta<0)){
    stop("Parameters not allowed")}
  random_genpois <- numeric(n)
  for (i in 1:n) {
    temp_random_genpois <- 0
    random_number <- runif(1)
    kum <- dgenp(0, lambda = lambda, theta = theta)
    while (random_number > kum) {
      temp_random_genpois <- temp_random_genpois + 1
      kum <- kum + dgenp(temp_random_genpois,
                         lambda = lambda, theta = theta)
    }
    random_genpois[i] <- temp_random_genpois
  }
  return(random_genpois)
}
