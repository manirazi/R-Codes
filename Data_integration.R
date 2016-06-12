closeAllConnections()
rm(list=ls())     #resetting all values
cat("\014")  #cleaning the screen

# example command to run the code in the terminal:
#Rscript --vanilla Trapz_data2.R -t data_xf.csv -t data_xf2.csv 


trapz_data<-function(data){  #function: trapezoidal integration method
  M<-nrow(data)-1     #number of grid points M
  SUM<-0
  for (K in 1:M){   #loop over subintervals
    SUM<-SUM+0.5*(data$x[K+1]-data$x[K])*(data$f[K+1]+data$f[K])
  }
  print("---trapezoidal integration method---")
  print(paste0("The approximate value of integral over the provided data set is: ",SUM))
  return(SUM)
}

# this function performed a numerical integration over a 
#randomly disributed data set
ScatteredDataIntegration<-function(data){
  M<-nrow(data)     #number of grid points M  
  ia<-1             #Points index starts from 1
  ib<-M             # ends at M
  n<-M
  if (n<4 | ia>ib){
    stop('the integration method failed', call=FALSE)
    if (ia>ib){
      stop('Usage: do not use an empty data file', call=FALSE)
    }
    if (n<4){
      stop('Usage: the number of data point for the filenames 
            must be greater than 4', call=FALSE)
    }
  } else if (ia==ib){
    result<-0.0
    error<-0.0
  } else{
    INT<-0.0
    error<-0.0
    s<-0.0
    c<-0.0
    r4<-0.0
    if (ia==n-1 && M>4){
      j<-n-2
    } else if (ia>2){
      j<-ia
    } else {
      j<-3
    }
    if (ib==1 && n>4){
      k<-4
    } else if (n>ib+2){
      k<-ib+1
    } else {
      k<-n-1
    }
    for (i in j:k){
      if (i==j){
        h2<-data$x[j-1]-data$x[j-2]
        d3<-(data$f[j-1]-data$f[j-2])/h2
        h3<-data$x[j]-data$x[j-1]
        d1<-(data$f[j]-data$f[j-1])/h3
        h1<-h2+h3
        d2<-(d1-d3)/h1
        h4<-data$x[j+1]-data$x[j]
        r1<-(data$f[j+1]-data$f[j])/h4
        r2<-(r1-d1)/(h4+h3)
        h1<-h1+h4
        r3<-(r2-d2)/h1
        if (ia==1){
          INT<-h2*(data$f[1]+h2*(0.5*d3-h2*((d2/6.0)-(h2+2.0*h3)*(r3/12.0))))
          s<- -(h2*h2*h2)*(h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
        }
      }
      if (i!=j){
        h4<-data$x[i+1]-data$x[i]
        r1<-(data$f[i+1]-data$f[i])/h4
        r4<-h4+h3
        r2<-(r1-d1)/r4
        r4<-r4+h2
        r3<-(r2-d2)/r4
        r4<-r4+h1
        r4<-(r3-d3)/r4
      }
      if (i<=ib && i>ia){
        INT<-INT+h3*((data$f[i]+data$f[i-1])*0.5-
                       h3*h3*(d2+r2+(h2-h4)*r3)/12.0)
        c<-(h3*h3*h3)*(2.0*h3*h3+5.0*(h3*(h4+h2)+2.0*h4*h2))/120.0
        error<-error+(c+s)*r4
        if (i==j){
          s<-s+2.0*c
        } else {
          s<-c
        }
      } else {
        error<-error+r4*s
      }
      if (i==k){
        if (ib==n){
          INT<-INT+h4*(data$f[n]-h4*(r1*0.5+
                                           h4*(r2/6.0+(2.0*h3+h4)*
                                                 r3/12.0)))
          error<-error-(h4*h4*h4)*r4*(h4*(3.0*h4+5.0*h2)+
                                        10.0*h3*(h2+h3+h4))/60.0
        }
        if (ib>=n-1){
          error=error+s*r4
        }
      } else {
        h1<-h2
        h2<-h3
        h3<-h4
        d1<-r1
        d2<-r2
        d3<-r3
      }
    }
  }
  result<-INT+error
  print("---highly accurated integration method on scattered data---")
  print(paste0('The approximate value of integral over the provided data set is: ',result))
  return(result)
}







library("optparse")
option_list=list(
  make_option (c("-t","--trapz"),type="character",default=NULL,
               help="trapezoidal integration method",
               metavar="character"),
  make_option (c("-s","--scatered"),type="character",default=NULL,
               help="highly accurated integration method on scattered data",metavar="character")
);

opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);
Number_of_files<-0
if (is.null(opt$t) && is.null(opt$s)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (integration method)", call=FALSE)
} else if (!is.null(opt$t)){
  args=commandArgs(trailingOnly = TRUE)
  if (length(args)==1){
    stop("at least one argument must be supplied (input 
         file) ", call=FALSE)
  }else if (length(args)>=2) {
    M<-length(args)
    filesToProcess<-args[2:M]
    K<-1
    print("++++NUMERICAL INTEGRATION ON DATA FILE++++")
    for (y in filesToProcess){
      if (y!="-t" && y!="-s"){
        Number_of_files<-Number_of_files+1
        print(paste0(Number_of_files,"--For data file ",y," :"))
        my_data<-read.csv(y,header = FALSE,col.names = c("x","f"))
        my_data<-my_data[complete.cases(my_data),]
        my_data<-my_data[with(my_data,order(x)),]
        rownames(my_data)<-1:nrow(my_data)
        if (args[K]=="-t"){
          A<-trapz_data(my_data)
        } else {
          A<-ScatteredDataIntegration(my_data)
        }
      }
      K<-K+1
    }
  }
} else {
  args=commandArgs(trailingOnly = TRUE)
  if (length(args)==1){
    stop("at least one argument must be supplied (input 
         file) ", call=FALSE)
  }else if (length(args)>=2) {
    M<-length(args)
    filesToProcess<-args[2:M]
    K<-1
    print("++++NUMERICAL INTEGRATION ON DATA FILE++++")
    for (y in filesToProcess){
      if (y!="-t" && y!="-s"){
        Number_of_files<-Number_of_files+1
        print(paste0(Number_of_files,"--For data file ",y," :"))
        my_data<-read.csv(y,header = FALSE,col.names = c("x","f"))
        my_data<-my_data[complete.cases(my_data),]
        my_data<-my_data[with(my_data,order(x)),]
        rownames(my_data)<-1:nrow(my_data)
        if (args[K]=="-t"){
          A<-trapz_data(my_data)
        } else {
          A<-ScatteredDataIntegration(my_data)
        }
      }
      K<-K+1
    }
  }
}


