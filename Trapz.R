func<- function(x) {  #defining a function to integrate
  f<- x*x
  return(f)
}
      #this program perform a one dimensional
      #numerical integration over [A,B] on a function
trapz<-function(A,B,M){
  H<-as.double(abs((B-A)/M))  #computing the spatial increment
  SUM<-0.0
  for (K in 1:M-1){     #loop over subintervals
    X<-A+as.double(K)*H
    SUM<-SUM+func(X)
  }
  SUM<-H*(func(A)+func(B)+2.0*SUM)/2.0  #computing the trapezoidal sum
  print(paste0('The approximate value of the integral of your desired function on the interval ',
               A,' ==> ',B,' using ',M,' subintervals is ',SUM))
  return(SUM)
}

A<-readline("Please enter the location of the 1D domain left bound A: ")
B<-readline("Please enter the location of the 1D domain right bound B: ")
M<-readline("Please enter your desired number of subintervals M: ")
A<-as.numeric(unlist(strsplit(A,",")))
B<-as.numeric(unlist(strsplit(B,",")))
M<-as.numeric(unlist(strsplit(M,",")))
trapz(A,B,M)
