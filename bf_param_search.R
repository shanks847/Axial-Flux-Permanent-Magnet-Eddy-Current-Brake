# Author :  Shankar Ramharack
# SID   :   816019096 
# Title :   Brute Force Paramater Search for APM-ECBS Design
#
# All functions use SI inputs hence they must be converted
# apriori. 
# ======================================================

# Define known variables
mu_0 <- (pi*4)*(1e-7) #permeability of air
ecp_ec <- (5.998e-7) #electrical conductivity of ecp
p = 4 #number of pole pairs
N <- 1000 #Upper limit of summation for fourier terms
K <- N
Br <- 1.25 #Magnetic remenance
pp <- 90 #pole-pitch
pa <- 40 #pole-arc
alpha <- pa/pp #pole-pitch to pole-arc ratio
a <- 5e-3 #APM mount thickness
e <- e  #ecp mount thickness
c <- 1e-3 #air gap thickness
d <- 2e-3 #ecp thickness

#Develop helper functions
solve_gamma_nk <- function(n, R_3,k,tau,ec,permeability_of_fs,ss,R_m) 
{
    return(sqrt(   ( ((n*pi)/(R_3))^2 + ((k*pi)/(tau))^2 + (1i *ec*permeability_of_fs*ss*R_m*(k*pi)/(tau)) )   ))
}

solve_a_nk <- function(n,R_3,k,tau)
{
    return(sqrt (((n*pi)/(R_3))^2 + ((k*pi)/(tau))^2 ))   
}

solve_r_underline <- function(a_nk,gamma_nk,b,c,d)
{
    numerator <- -(cosh(a_nk *c)*sinh(gamma_nk*d) + ((a_nk)/(gamma_nk))*sinh(a_nk *c)*cosh(gamma_nk*d) )
    denominator <- (cosh(a_nk *(b+c))*sinh(gamma_nk*d) + ((a_nk)/(gamma_nk))*sinh(a_nk *(b+c))*cosh(gamma_nk*d) )
    return(numerator/denominator)
        
}

solve_R_m <- function(R_1, R_2)
{
    return ((R_1 + R_2)/2)
}

solve_M_nk <- function(mag_rem, permeability_of_fs,n,k,alpha, R_1,R_2,R_3)
{
    t1 <- (16*mag_rem)/((pi^2) * permeability_of_fs*n*k) #first term
    t2 <- sin(k*alpha*pi/2) #second term
    t3 <- sin(n*(pi/2)*((R_2 - R_1)/R_3)) #third term
    return(t1*t2*t3)
}

#Torque equation
solve_em_torque <- function(em_vec,radii_vec,geo_vec,N,K,ss)
{
    #  To make input less cumbersome it takes 3 vectors as well as indiv args
    #  To ensure correct results please ensure the vecs sequence is maintained
    #  when passing values
    #  
    #  em_vec <- c(permeability of fs, mag_rem, ec)
    #  radii_vec <- c(R_1, R_2, R_3)
    #  geo_vec <- c(b,c,d,p,alpha)
    #
    
    R_m <- solve_R_m(radii_vec[1],radii_vec[2])
    tau <- 0.5*pi*R_m
    t1 <- 0.5*em_vec[1]*(geo_vec[4]^2)*tau*radii_vec[3]
    t2 <- 0
    
    for (n in 1:N) 
    {
        for (k in 1:K)
        {
            gamma_nk <- solve_gamma_nk(n,radii_vec[3],k,tau,em_vec[3],em_vec[1],ss,R_m)
            M_nk <- solve_M_nk(em_vec[2],em_vec[1],n,k,geo_vec[5],radii_vec[1],radii_vec[2],radii_vec[3])
            a_nk <- solve_a_nk(n,radii_vec[3],k,tau)
            r_underline <- solve_r_underline(a_nk,gamma_nk,geo_vec[1],geo_vec[2],geo_vec[3])
            t2 <- t2 + (1i *k*((M_nk^2)/(a_nk)) * r_underline*sinh(a_nk *geo_vec[1]) )
        }
    }
    
    return(t1*Re(t2))
}

############ TESTING #####################
t_em_vec <- c(mu_0,Br,ecp_ec)
t_radii_vec <- c(0.12,0.1,0.15)
t_geo_vec <- c(0.02,c,d,p,alpha)
t_Te <- solve_em_torque(t_em_vec,t_radii_vec,t_geo_vec,N,K,2000)



#TODO Modify the code to plot Torque results over parameter space
library(plotly)


DF <- data.frame(a = 1:3, b = letters[10:12],
                 c = seq(as.Date("2004-01-01"), by = "week", length.out = 3),
                 stringsAsFactors = TRUE)
data.matrix(DF[1:2])
#data.matrix(DF)
fig <- plot_ly(z = ~volcano)
fig <- fig %>% add_surface()
fig
