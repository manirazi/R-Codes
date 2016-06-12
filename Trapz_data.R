closeAllConnections()
rm(list=ls())     #resetting all values
cat("\014")  #cleaning the screen

#This program works on the R terminal and 
#computes the numerical integration over a data set

trapz_data<-function(data){
  M<-nrow(data)-1     #number of grid points M
  SUM<-0
  for (K in 1:M){   #loop over subintervals
    SUM<-SUM+0.5*(data$x[K+1]-data$x[K])*(data$f[K+1]+data$f[K])
  }
  print(paste0("The approximate value of integral over the provided data set is: ",SUM))
  return(SUM)
}


filesToProcess<-dir(pattern = "*\\.csv")  #find .csv files 
                                           #in the working directory
K<-1
print("++++NUMERICAL INTEGRATION ON DATA FILE++++")
for (y in filesToProcess){
  print(paste0(K,"--For data file ",y," :"))
  # reafing data files with assagining column names
  my_data<-read.csv(y,header = FALSE,col.names = c("x","f"))
  #removing rows with NaNs
  my_data<-my_data[complete.cases(my_data),]
  #sorting the data
  my_data<-my_data[with(my_data,order(x)),]
  #reordering column indices 
  rownames(my_data)<-1:nrow(my_data)
  A<-trapz_data(my_data)
  K<-K+1   #file number
}

