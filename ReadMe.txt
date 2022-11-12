
PROJECT TITLE:
Outlier Detection in Non-Parametric Profile Monitoring




Presence of outliers in the data highly affect the performance of the various models. In our project we replicated a paper which deals with this problem of outliers using a non-parametric approach. We replicated the two-stage algorithm in R software and then simulated data from the two settings – Auto Regressive and Moving Average  in order to test the performance of our algorithm. 
We have used a non-parametric outlier detection measure based on Kernel function. In this Kernel function we needed an optimum bandwidth and we tried to find it using 2 methods – optimized function and by plotting the errors. Then used Least Trimmed Kernel Distance algorithm to find the detection measure.
This algorithm uses an iteration process to find the ultimate clean subset and then using hypothesis testing, it tries to eliminate the outlying profiles in data.
The profiles using LTKD procedure may not be reliable so we used a refinement scheme to enhance the detection efficiency.

The code performs the above mentioned method in sequence and provide the RLTKD value in the end which suggests which should be the outlying and non-outlying profile. 







