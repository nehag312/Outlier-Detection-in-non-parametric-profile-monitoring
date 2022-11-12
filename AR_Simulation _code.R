
#*******************************SIMULATION STUDIES USING AUTO-REGRESSIVE DATA*************************************************


#Function to generate covariates using Auto-Regressive method

covariate_ar<-function(m,n)
{
  # mean vector
  mean1 = c(0,0,0,0)
  # covariance matrix
  sigma = matrix(0,nrow=4,ncol=4)
  for(i in seq(1:4)){
    for(j in seq(1:4)){
         sigma[i,j] = (0.5)**(abs(i-j))
                       }
                     }
  xij<-vector("list",m)
  for(i in 1:m){
    xij[[i]]=matrix(0,nrow=n,ncol=4)
    for (j in 1:n){
      xij[[i]][j,]=mvrnorm(1,mean1,sigma)
                  }
                }
  return(xij)
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

# We are using m = 50 ( number of profiles) and n=100 (number of response vector for each profile)

### Covariates are of dimension = 4

m1<-10    # Globally defined
mo<-40
m=m1+mo
n<-200
ano=1.1
ao=2


# Simulation of covariates in case of AR structure

set.seed(29)
data_x_ar = covariate_ar(m,n)
data_x_ar

# Generation of response matrix in case of AR structure 

# Matrix of gx in case of AR

mat_gx_ar<-matrix(0,m,n)
for(i in 1:mo){
  for(j in 1:n){
    mat_gx_ar[i,j]=gaxij(ao,data_x_ar[[i]][j,])
  }
}
for(i in (mo+1):m){
  for(j in 1:n){
    mat_gx_ar[i,j]=gaxij(ano,data_x_ar[[i]][j,])
  }
}

# Matrix of epsilon

eij_ar<-matrix(rnorm(m*n,0,1),m,n)

# Response matrix in case of AR structure

data_y_ar<-mat_gx_ar+eij_ar
data_y_ar

# Plotting of profiles

plot(1:100,mat_gx_ar[1,1:100],"l",lwd=1.5,ylim=c(-10,300),col="green",main="Outlying and Non-Outlying Profile ",ylab="y",xlab="X")
for(i in 2:m1){
  lines(1:100,data_y_ar[i,1:100],"l",col="green")
              }
for(i in (m1+1):m){
  lines(1:100,data_y_ar[i,1:100],"l")
                  }


# We first need to divide the data into two parts 

## Part 1. To obtain the clean subset of profiles - used for estimating our gij_hat
## Part 2. To do hypothesis testing - to determine which all profiles are outlying profiles 


#Just for the sake of simplicity, let's divide data_y_ar and data_x_ar as follows :

data_y_ar1<-data_y_ar[,1:100]
data_y_ar2<-data_y_ar[,101:200]  

data_x_ar1<-vector("list",m);data_x_ar2<-vector("list",m)
for(i in 1:m){
  data_x_ar1[[i]]=data_x_ar[[i]][1:100,]
  data_x_ar2[[i]]=data_x_ar[[i]][101:200,]
             }


# In order to find Detection measure, we define kernel function which takes w as bandwidth (will find later)

# kernel function (Gaussian)

kernel<-function(w,i,k,j){
  const=((1/(2*pi))^(2))*(1/sqrt(det(diag(w,4))))
  x=as.matrix((data_x_ar1[[i]][k,]-data_x_ar1[[i]][j,]))
  k=const*exp(-as.numeric(t(x)%*%x)/(2*w))
  return(k)
                        }


# OPTIMAL BANDWIDTH FOR THE ITH PROFILE

# Estimation of gij( ith regression function by eliminating the jth variable )


gi_j<-function(w,i){
  mat_gi_j<-c()
  for(j in 1:(n/2)){
    num=0
    den=0
    for(k in 1:(n/2)){
      if(k!=j){
        num=num+data_y_ar1[1,k]*kernel(w,i,k,j)
        den=den+kernel(w,i,k,j)
      }
    }
    mat_gi_j[j]=num/den}
  return(mat_gi_j) }

#*******************************************************************
# Using optimize function 

theta=1
ei_fun<-function(par=w){
gi_j_hat=gi_j(par,theta)
return(mean((data_y_ar1[theta,]-gi_j_hat)^2))

wih<-c()
for(i in 1:m){
   wih[i]<-optimize(ei_fun,interval=c(0.001,2),maximum=F)$minimum
   theta=theta+1
             }
median(wih)   ##we see that the optimal bandwidth is occurring in the extreme. 

#********************************************************************

# So, we have used Silvermann's Optimal Bandwith

opt.bandwidth<-function(sample){
h<-c()
for(i in 1:ncol(sample)){
h[i]<-1.06*sqrt(var(sample[,i]))*length(sample[,i])^(-1/5)}
return(h)
                              }
wi<-c()
for(i in 1:m){
  wi[i]<-mean(opt.bandwidth(data_x_ar[[i]]))
             }
wi
opt_h<-min(wi)
#*********************************************************************

# Then we tried using plottoing the data and found opt_h = 2.5(approx)

opt_h=2.5

# Initial subset of H={1,......,m} to start the iteration (H_int)

h_vec = sample(1:50,2,replace=TRUE)
h_vec

# OBTAINING A CLEAN PROFILE SET


# Finding estimate of gij_hat using h_vec (H_int) values as indexes

gij_h_hat<-function(i,j,h_vec){
  num<-0
  den<-0
  for(h in h_vec){
    num<-num+data_y_ar1[i,h]*kernel(opt_h,i,h,j)
  }                                                           
  for(h in h_vec){
    den<-den+kernel(opt_h,i,h,j)
  }
  gij_hat = num/den
  return(gij_hat)
}

Finding estimate of gij_hat using h_vec (H_int) values as indexes

# Function for estimating Bias values 

cij_h_hat<-function(i,j,h_vec){
  c<-data_y_ar1[i,j]-gij_h_hat(i,j,h_vec)
  return(c)
}

# Function for finding Outlier Detection Measure for ith profile

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

detection_measure = function(h_vec) 
{
  din = c()
  for(i in 1:m) {
    din[i]= din(i,h_vec)
  }
  return(din)
}

# Finding optimal h indexes ( After sorting squares of detection measures )

h_index<-function(det_measure){
  srt=sort(det_measure^2)
  ind = c()
  for(i in 1:(floor(length(det_measure)/2)+1)){
    ind[i]=match(srt[i],det_measure^2)
  }
  return(ind)
}

# Function to check the convergence to obtain the ultimate clean subset

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

#******************************************************************************************************
#  First Iteration 

det_measure_1 = detection_measure(h_vec)
det_measure_1

h1= h_index(det_measure_1)
h1

#  Second Iteration 

det_measure_2 = detection_measure(h1)
det_measure_2

h2 = h_index(det_measure_2)
h2

# Convergence check

convergence_check(det_measure_1,det_measure_2,h1,h2)

det_measure_3 = detection_measure(h2)
det_measure_3

h3 = h_index(det_measure_3)
h3

convergence_check(det_measure_2,det_measure_3,h2,h3)

### We do it until we get the result "CONVERGED" from the convergence function
### Suppose we converged at h3, so we set h_ltkd = h3 (clean subset)

h_ltkd<-h3

### HYPOTHESIS TESTING USING PART - 2 of the data

# We redefined all the functions for part2 of the data

# kernel function (Gaussian) for part2 data

kernel_part2<-function(w,i,k,j){
  const=((1/(2*pi))^(2))*(1/sqrt(det(diag(w,4))))
  x=as.matrix((data_x_ar2[[i]][k,]-data_x_ar2[[i]][j,]))
  k=const*exp(-as.numeric(t(x)%*%x)/(2*w))
  return(k)
}

# gij_hat using h_vec valeues as indexes for part 2 data

gij_h_hat_part2<-function(i,j,h_vec){
  num<-0
  den<-0
  for(h in h_vec){
    num<-num+data_y_ar2[i,h]*kernel_part2(opt_h,i,h,j)
  }                                                           
  for(h in h_vec){
    den<-den+kernel_part2(opt_h,i,h,j)
  }
  gij_hat = num/den
  return(gij_hat)
}

# Biases for part2 data

cij_h_hat_part2<-function(i,j,h_vec){
  c<-data_y_ar2[i,j]-gij_h_hat_part2(i,j,h_vec)
  return(c)
}

# Outlier detection measure for ith profile for part 2 data

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


# Now, we find the Outlier Detection Measure for all profiles using h_ltkd (Clean subset)

det_measure_ltkd = detection_measure_part2(h_ltkd)
det_measure_ltkd


#******************Using Proposition 2.2 from the paper, we find find the asymptotic distribution of[ n*(opt_h)^2 * din(h_ltkd)]*****************

# Calculating estimate of its covariance matrix

sigma_i_est <-function(i,h_vec){
  ksum=0
  for(k in 1:(n/2)){
    for(l in 1:(n/2)){
      if(l!=k){
        ksum=ksum+((kernel_part2(opt_h,i,k,l))**2)*((cij_h_hat_part2(i,l,h_vec))**2)*((cij_h_hat_part2(i,k,h_vec))**2)
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

# Now we remove Outlier profiles based on threshold

# we take alpha = 0.5 

for(i in 1:m) {
  if(stand_din_h_ltkd[i] > 1.64 | stand_din_h_ltkd[i] < -1.64) {
    cat(" The ",i,"-th profile in our generated dataset is deemed as an outlier","\n")
              }
else{ cat(" The ",i,"-th profile in our generated dataset is not an outlier","\n") }
}


##****************************REFINEMENT STEP*************************



# Compliment of outlying profiles set

t = seq(1:m) # Set of all the profiles

# Suppose we obtained h, a set containing the profiles that came out to be outlier in the previous step
# The find a compliment of that i.e. t\h = h_c
  

# Finding Detection Measure using h_c

det_measure_ref = detection_measure_part2(h_c)
det_measure_ref

# Finding Standardized Detection Measure using h_c

stand_din_h_ref = T_h_part2(h_c)
stand_din_h_ref

# we take alpha = 0.5 

for(i in 1:m) {
  if(stand_din_h_ref[i] > 1.64 | stand_din_h_ref[i] < -1.64) {
    cat(" The ",i,"-th profile in our generated dataset is deemed as an outlier","\n")
                                                             }
  else{ cat(" The ",i,"-th profile in our generated dataset is not an outlier","\n") }
              }

# Results from above step give us a SET OF OUTLYING PROFILES

############################################# END OF CODE #################################################
 
