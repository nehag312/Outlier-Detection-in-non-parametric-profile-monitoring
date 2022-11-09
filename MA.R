library(MASS)

# Function to generate covariates using Moving Average method
#it will be of four dimension
#NOTE- this function doesn't need any specific input from the user

covariate<-function(x){
  #generate two neta rv from U(0,1) independently
  neta<-c(runif(2,0,1))
  
  #generate 5 zijk from N|(0,1)
  zijk<-c(rnorm(5,0,1))
  
  xijk<-c()
  for(i in 1:4){
    xijk[i]=(neta[1]*zijk[i]+neta[2]*zijk[i+1])/sqrt(sum(neta^2))
  }
  return(xijk )
}

# Borel measurable function g(.), which determines the relationship between covariates and response

gaxij<-function(a,covariate){
  covariate<-as.vector(covariate)
  gsum=a*covariate[1]+(2*a*covariate[2]-1)^2+
    sin(2*pi*a*covariate[3])/(2-sin(2*pi*a*covariate[3]))+
    0.1*sin(2*pi*a*covariate[4])+0.2*cos(2*pi*a*covariate[4])+
    0.3*sin((2*pi*a*covariate[4])^2)+0.4*cos((2*pi*a*covariate[4])^3)+
    0.5*sin((2*pi*a*covariate[4])^3)
  return(gsum)
}

## We are using m = 50 (number of profiles) AND n=100 (number of response vector for each profile)*

# Globally defined
m1<-45
mo<-5
m=m1+mo
n<-100
ano=1.1
ao=2


##Simulation of covariates

#Covariates generated in case of MA structure

set.seed(29)

#from outlying profiles
data_x_ma<-vector("list",m)
for(i in 1:m){
  data_x_ma[[i]]=matrix(0,nrow=n,ncol=4)
  for (j in 1:n){
    data_x_ma[[i]][j,]=covariate(j)  # used function covariate
  }
}

##Generation of response matrix

#Response matrix in case of MA structure 

# used global data_x_ma

# Matrix of gx in case of MA

mat_gx_ma<-matrix(0,m,n)
for(i in 1:mo){
  for(j in 1:n){
    mat_gx_ma[i,j]=gaxij(ao,data_x_ma[[i]][j,])
  }
}
for(i in (mo+1):m){
  for(j in 1:n){
    mat_gx_ma[i,j]=gaxij(ano,data_x_ma[[i]][j,])
  }
}

# Matrix of epsilon
set.seed(29)
eij_ma<-matrix(rnorm(m*n,0,1),m,n)

# Response matrix in case of MA structure

data_y_ma<-mat_gx_ma+eij_ma


data_y_ma

##plotting the profiles
#coloring the outlying profiles in red

par(mfrow=c(1,1))
plot(1:100,mat_gx_ma[1,1:100],"l",lwd=1.5,ylim=c(-10,300),col="red",main="Outlying and Non-Outlying Profile ",ylab="y",xlab="X")
for(i in 2:5){
  lines(1:100,data_y_ma[i,1:100],"l",col="red")
}
for(i in 6:50){
  lines(1:100,data_y_ma[i,1:100],"l")
}

##We first need to divide the data into two parts

##1.   To obtain the clean subset of profiles - used for estimating our gij_hat
##2.   To do hypothesis testing - to determine which all profiles are outlying profiles 

#Just for the sake of simplicity, let's divide data_y_ma and data_x_ma as follows
data_y_ma1<-data_y_ma[,1:50]
data_y_ma2<-data_y_ma[,51:100]  #will be used later

data_x_ma1<-vector("list",m);data_x_ma2<-vector("list",m)
for(i in 1:m){
  data_x_ma1[[i]]=data_x_ma[[i]][1:50,]
  data_x_ma2[[i]]=data_x_ma[[i]][51:100,]
}

##In order to find Detection measure, we define kernel function which takes w as bandwidth (will find later)
#here we have used the Gaussian kernel

kernel<-function(w,i,k,j){
  const=((1/(2*pi))^(2))*(1/sqrt(det(diag(w,4))))
  x=as.matrix((data_x_ma1[[i]][k,]-data_x_ma1[[i]][j,]))
  k=const*exp(-as.numeric(t(x)%*%x)/(2*w))
  return(k)
}

##*Finding optimal bandwidth
#gi_j is used for cross validation 
#It is Nadaraya-Watson Kernel estimator of the ith regression function by omitting the jth variable
gi_j<-function(w,i){
  mat_gi_j<-c()
  for(j in 1:(n/2)){
    num=0
    den=0
    for(k in 1:(n/2)){
      if(k!=j){
        num=num+data_y_ma1[1,k]*kernel(w,i,k,j)
        den=den+kernel(w,i,k,j)
      }
    }
    mat_gi_j[j]=num/den}
  return(mat_gi_j) }

#OPTIMAL BANDWIDTH FOR THE ITH PROFILE
shagun=1
ei_fun<-function(par=w){
  gi_j_hat=gi_j(par,shagun)
  return(mean((data_y_ma1[shagun,]-gi_j_hat)^2))
}
wih<-c()
for(i in 1:m){
  wih[i]<-optimize(ei_fun,interval=c(0.001,2),maximum=F)$minimum
  shagun=shagun+1
}
median(wih)   ##we see that the optimal bandwidth is occurring in the extreme.

###plotting ei_fun against w
par(mfrow=c(2,5))
ei_val=c()
w<-seq(0.01,7,length.out=100)
shagun=1
while(shagun<21){
  for(i in 1:100){ ei_val[i]=ei_fun(w[i])}
  plot(w,ei_val)
  shagun=shagun+1}

#we have also tried to use Silvermann's Optimal Bandwidth but didn't get good results from that
opt.bandwidth<-function(sample){
  h<-c()
  for(i in 1:ncol(sample)){
    h[i]<-1.06*sqrt(var(sample[,i]))*length(sample[,i])^(-1/5)}
  return(h)
}
wi<-c()
for(i in 1:m){
  wi[i]<-mean(opt.bandwidth(data_x_ma[[i]]))
}

# opt_h<-median(wi)  #in case of Silvermanns
# opt_h<-2.5         #in case of minimising the squared error 

#taking only two profiles
#h_vec = sample(1:m,2,replace=TRUE)
h_vec = c(14,27)

## gij_hat using h_vec values as indexes
gij_h_hat<-function(i,j,h_vec){
  num<-0
  den<-0
  for(h in h_vec){
    num<-num+data_y_ma1[i,h]*kernel(opt_h,i,h,j)
  }                                                           
  for(h in h_vec){
    den<-den+kernel(opt_h,i,h,j)
  }
  gij_hat = num/den
  return(gij_hat)
}

##the estimated bias value
cij_h_hat<-function(i,j,h_vec){
  c<-data_y_ma1[i,j]-gij_h_hat(i,j,h_vec)
  return(c)
}

# outlier detection measure for ith profile

din<-function(i,h_vec){
  ksum=0
  for(k in 1:(n/2)){
    for(l in 1:(n/2)){
      if(l!=k){
        ksum=ksum+kernel(opt_h,i,k,l)*cij_h_hat(i,l,h_vec)*cij_h_hat(i,k,h_vec)
      }
    }
  }
  return(ksum/((n/2)*(n/2-1)))
}


# Outlier detection measure for all profiles


#we constructed the detection meaure for all profile using the above function in it
detection_measure = function(h_vec) 
{
  din = c()
  for(i in 1:m) {
    din[i]= din(i,h_vec)
  }
  return(din)
}

# finding optimal h indexes(which most probably do not contain indexes of outlying profile)
h_index<-function(det_measure){
  srt=sort(det_measure^2)
  ind = c()
  for(i in 1:(floor(length(det_measure)/2)+1)){
    ind[i]=match(srt[i],det_measure^2)
  }
  return(ind)
}

##To test if the Din square converges or not
convergence_check = function(det_measure1,det_measure2,h1,h2)
{
  sum1=0
  sum2=0
  for(i in seq_along(h1)) {
    sum1 = sum1+det_measure1[h1[i]]**2
  }
  for(j in seq_along(h2)) {
    sum2 = sum2+det_measure2[h2[j]]**2
  }
  if(abs(sum1-sum2)<0.001) { print("converged") }
  else { print("Not converged yet, move to next iteration") }
}

##First Iteration 

det_measure_1 = detection_measure(h_vec)
det_measure_1

h1= h_index(det_measure_1)
h1

det_measure_2 = detection_measure(h1)
det_measure_2

h2 = h_index(det_measure_2)
h2

convergence_check(det_measure_1,det_measure_2,h1,h2)

det_measure_3 = detection_measure(h2)
det_measure_3

h3 = h_index(det_measure_3)
h3

convergence_check(det_measure_2,det_measure_3,h2,h3)
#this convereger

#ultimate clean subset
h_ltkd<-h3

# Suppose the clean subset is obtained , h_ltkd

# kernel function (Gaussian) for part2 data
kernel_part2<-function(w,i,k,j){
  const=((1/(2*pi))^(2))*(1/sqrt(det(diag(w,4))))
  x=as.matrix((data_x_ma2[[i]][k,]-data_x_ma2[[i]][j,]))
  k=const*exp(-as.numeric(t(x)%*%x)/(2*w))
  return(k)
}

## gij_hat using h_vec valeues as indexes for part 2 data
gij_h_hat_part2<-function(i,j,h_vec){
  num<-0
  den<-0
  for(h in h_vec){
    num<-num+data_y_ma2[i,h]*kernel_part2(opt_h,i,h,j)
  }                                                           
  for(h in h_vec){
    den<-den+kernel_part2(opt_h,i,h,j)
  }
  gij_hat = num/den
  return(gij_hat)
}

# biases for part2 data

cij_h_hat_part2<-function(i,j,h_vec){
  c<-data_y_ma2[i,j]-gij_h_hat_part2(i,j,h_vec)
  return(c)
}

# outlier detection measure for ith profile for part 2 data

din_part2<-function(i,h_vec){
  ksum=0
  for(k in 1:(n/2)){
    for(l in 1:(n/2)){
      if(l!=k){
        ksum=ksum+kernel_part2(opt_h,i,k,l)*cij_h_hat_part2(i,l,h_vec)*cij_h_hat_part2(i,k,h_vec)
      }
    }
  }
  return(ksum/((n/2)*(n/2-1)))
}

# Outlier detection measure for all profiles of part 2 data

detection_measure_part2 = function(h_vec) 
{
  din = c()
  for(i in 1:m) {
    din[i]= din_part2(i,h_vec)
  }
  return(din)
}



# Now its time to use part 2 of data  for hypothesis testing

# Outlier detection measure for all profiles using h_ltkd
det_measure_ltkd = detection_measure_part2(h_ltkd)
det_measure_ltkd

# Using Proposition 2.2 , we find find the asymptotic distribution of n*(opt_h)^2 * din(h_ltkd)

# Calculating estimate of its covariance matrix

sigma_i_est <-function(i,h_vec){
  ksum=0
  for(k in 1:(n/2)){
    for(l in 1:(n/2)){
      if(l!=k){
        ksum=ksum+(kernel_part2(opt_h,i,k,l))*(cij_h_hat_part2(i,l,h_vec)**2)*(cij_h_hat_part2(i,k,h_vec)**2)
      }
    }
  }
  return(2*ksum/((n/2)*(n/2-1)))
}

# Finding the standardized version of Outlier Detection Measure Din(h_ltkd)

n1=n/2
T_h_part2 = function(h_ltkd) {
  ti_h_ltkd = c()
  for(i in 1:m) {
    const = sqrt((n1-1)/n1)*n1*(opt_h)**2
    ti_h_ltkd[i]= const*((din_part2(i,h_ltkd))/sqrt(sigma_i_est(i,h_ltkd)))
  }
  return(ti_h_ltkd)           
}

stand_din_h_ltkd = T_h_part2(h_ltkd)
stand_din_h_ltkd

# Removing Outlier profiles based on threshold

# we take alpha = 0.5 

for(i in 1:m) {
  if(stand_din_h_ltkd[i] > 1.64 | stand_din_h_ltkd[i] < -1.64) {
    cat(" The ",i,"-th profile in our generated dataset is deemed as an outlier","\n")
  } else{print("No outlying profile")
  }
}

# *Refinement Step*
# Compliment of outlying profiles set

s = seq(1:m)
h = c(0)  #obtained in above step
h_c <-s #s[pmatch(s,h)  #compliment of h

#It is possible that still some outlying profile maybe present
#So we used the refined detection measure, for checking the same 
det_measure_ref = detection_measure_part2(h_c)
det_measure_ref


##standard version of the refined outlier detection measure
stand_din_h_ref = T_h_part2(h_c)
stand_din_h_ref
