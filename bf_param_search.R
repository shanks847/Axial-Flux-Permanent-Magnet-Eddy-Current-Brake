# Author :  Shankar Ramharack
# SID   :   816019096 
# Title :   Brute Force Paramater Search for APM-ECBS Design
#
# All functions use SI inputs hence they must be converted
# apriori. 
# ======================================================
library(progress)
library(dplyr)


# Define known variables
mu_0 <- (pi*4)*(1e-7) #permeability of air
ecp_ec <- 57.7e6 #electrical conductivity of ecp
p = 4 #number of pole pairs
N <- 400 #Upper limit of summation for fourier terms
K <- N
Br <- 1.25 #Magnetic remenance
pp <- 90 #pole-pitch
pa <- 40 #pole-arc
alpha <- pa/pp #pole-pitch to pole-arc ratio
a <- 5e-3 #APM mount thickness
#e <- e  #ecp mount thickness
c <- 1e-3 #air gap thickness
d <- 2e-3 #ecp thickness

#Develop helper functions
solve_gamma_nk <- function(n, R_3,k,tau,ec,permeability_of_fs,ss,R_m) 
{
    # n,k   := Fourier series index
    # R_3 := radius of conducting plate
    # tau := pi/pole_pairs_num * magnetic remenance
    # ec := electrical conductivity of conductor
    # ss := slip speed (RPM)
    # permeabiility_of_fs := permeability of free space
    # R_m := mean readius of magnet arrangement
    return(sqrt(   ( ((n*pi)/(R_3))^2 + ((k*pi)/(tau))^2 + (1i *ec*permeability_of_fs*ss*R_m*(k*pi)/(tau)) )   ))
}

solve_a_nk <- function(n,R_3,k,tau)
{
    # n,k   := Fourier series index
    # R_3 := readius of conducting plate
    # tau := pi/pole_pairs_num * magnetic remenance
    return(sqrt (((n*pi)/(R_3))^2 + ((k*pi)/(tau))^2 ))   
}

solve_r_underline <- function(a_nk,gamma_nk,b,c,d)
{
    # b := magnet thickness
    # c := airgap thickness
    # d := conducting plate thickness
    
    numerator <- -(cosh(a_nk *c)*sinh(gamma_nk*d) + ((a_nk)/(gamma_nk))*sinh(a_nk *c)*cosh(gamma_nk*d) )
    denominator <- (cosh(a_nk *(b+c))*sinh(gamma_nk*d) + ((a_nk)/(gamma_nk))*sinh(a_nk *(b+c))*cosh(gamma_nk*d) )
    return(numerator/denominator)
        
}

solve_R_m <- function(R_1, R_2)
{
    # R_1 := inner radius of magnet
    # R_2 := outer radius of magnet
    return ((R_1 + R_2)/2)
}

solve_M_nk <- function(mag_rem, permeability_of_fs,n,k,alpha, R_1,R_2,R_3)
{
    # R_1 := inner radius of magnet
    # R_2 := outer radius of magnet
    # R_3 := radius of conducting plate
    # permeabiility_of_fs := permeability of free space
    # mag_rem := magnetic remenance of magnet
    # alpha := pole-pitch to pole-arc ratio
    # n,k   := Fourier series index
    
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
    tau <- (1/geo_vec[3])*pi*R_m
    t1 <- 0.5*em_vec[1]*(geo_vec[4]^2)*tau*radii_vec[3]
    t2 <- 0
    #est_iterations = N*K
    #pb <- progress_bar$new(
    #    format = "  Walking parameter space [:bar] :percent | Elapsed :elapsed",
    #    total = est_iterations, clear = FALSE, width= 60)
    for (n in 1:N) 
    {
        for (k in 1:K)
        {
            gamma_nk <- solve_gamma_nk(n,radii_vec[3],k,tau,em_vec[3],em_vec[1],ss,R_m)
            M_nk <- solve_M_nk(em_vec[2],em_vec[1],n,k,geo_vec[5],radii_vec[1],radii_vec[2],radii_vec[3])
            a_nk <- solve_a_nk(n,radii_vec[3],k,tau)
            r_underline <- solve_r_underline(a_nk,gamma_nk,geo_vec[1],geo_vec[2],geo_vec[3])
            t2 <- t2 + (1i *k*((M_nk^2)/(a_nk)) * r_underline*sinh(a_nk *geo_vec[1]) )
            #pb$tick()
        }
    }
    #print("T2")
    #print(t2)
    
    return(t1*Re(t2))
}


walk_pspace <- function(em_vec,fixed_geo_vec,fixed_radii_vec,b_lims,wm_lims,N,K,ss)
{
    # em_vec = (permeability of free space, mag. remenance of conductor,
    #                               electrical conductivity of conductor)
    # fixed_radii_vec = (fixed R_2 and R_3 radii respectively)
    # fixed_geo_vec = (fixed c,d,p and alpha respectively)
    # b_lims = (lower bound , upper bound, step size) of magnet thickness
    # wm_lims = (lower bound , upper bound,step size ) of magnet radiual extrusion
    est_iterations =(length(seq(from=wm_lims[1], to = wm_lims[2], by = wm_lims[3])) * length(seq(from=b_lims[1], to = b_lims[2], by = b_lims[3])) )
    pb <- progress_bar$new(
        format = "  Walking parameter space [:bar] :percent | Elapsed :elapsed",
        total = est_iterations, clear = FALSE, width= 60)
    #pb$tick(0)
    
    R_2 = fixed_radii_vec[1]
    R_3 = fixed_radii_vec[2]
    R_1 = 0
    i_Te = 0
    completed_iter = 1
    
    tdf <- data.frame(matrix(ncol = 3, nrow = 0))
    #provide column names
    colnames(tdf) <- c('PM Axial Extrusion(w_m)', 'Magnet thickenss', 'Torque')
    #tdf[nrow(tdf) + 1,] <- c(10, 20, 30)
    
    for (wm in seq(from=wm_lims[1], to = wm_lims[2], by = wm_lims[3]))
    {
        R_1 = R_2 - wm
        i_radii_vec = c(R_1,R_2,R_3)
        #print(i_radii_vec)
        for (b in seq(from=b_lims[1], to = b_lims[2], by = b_lims[3]))
        {
            #print(paste0("[+]Iteration ", completed_iter,"  of ",est_iterations))
            pb$tick()
            #print(b)
            i_geo_vec = c(b,fixed_geo_vec[1],fixed_geo_vec[2],fixed_geo_vec[3],fixed_geo_vec[4])
            #print(i_geo_vec)
            i_Te = (1/23256) * solve_em_torque(em_vec,i_radii_vec,i_geo_vec,N,K,ss)#scaling
            
            #print(i_Te)
            #print( c(wm, b, i_Te))
            tdf[nrow(tdf) + 1,] <- c(wm, b, i_Te)
            #completed_iter = completed_iter + 1
        }
        
    }
    
    return(tdf)
    
    
}










############ TESTING ##################### 
t_em_vec <- c(mu_0,Br,ecp_ec)
t_radii_vec <- c(0.10,0.12,0.15)
t_geo_vec <- c(0.02,c,d,p,alpha)
#t_Te <- solve_em_torque(t_em_vec,t_radii_vec,t_geo_vec,N,K,2000)
#print(t_Te)

#constructing parameter domains for TORQUE by varying radius,magnet 

# for (x in seq(from=5, to = 100, by = 5))
# {
#     print(x)
# }

#t_wm_lims = c(0.02,0.11,0.01)
#t_b_lims = c(0.004,0.01,0.001)
t_wm_lims = c(0.02,0.11,0.01)
t_b_lims = c(0.005,0.04,0.005) #fontchagnester
# t_fixed_radii_vec = c(0.112,0.14) #Conducting disc is 140 mm, mount is 112 mm
# t_fixed_geo_vec = c(c,d,p,alpha)
# psp_df1k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,1000)
# psp_df2k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,2000)
# psp_df4k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,4000)
# psp_df6k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,6000)
# psp_df8k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,8000)
# 



#################### VERIFYING FONTCHAGNESTER #########################
v_em_vec <- c(mu_0,Br,ecp_ec)
v_radii_vec <- c(0.091,0.1315,0.15)
v_geo_vec <- c(0.00433,0.003,0.00421,10,0.71)
#vres <- solve_em_torque(v_em_vec,v_radii_vec,v_geo_vec,2000,2000,75)

#The result when first simulated was found to be scaled by 23259.14 Hence, to reduce error
#all subsequent calculations will be scaled by 1/232560. Futhermore, N = 1000 for calcs to reduce computation time
#Accuracy should be maintained

modified_psp_df1k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,1000)
modified_psp_df2k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,2000)
modified_psp_df4k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,4000)
modified_psp_df6k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,6000)
modified_psp_df8k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,N ,K,8000)


tmodified_psp_df1k <- walk_pspace(t_em_vec,t_fixed_geo_vec,t_fixed_radii_vec,t_b_lims,t_wm_lims,500 ,500,1000)
seq(from=t_b_lims[1], to = t_b_lims[2], by = t_b_lims[3])


######################## ANALYZING DATA ###########################333

psp1kres <- data.frame(modified_psp_df1k)
psp2kres <- data.frame(modified_psp_df2k)
psp4kres <- data.frame(modified_psp_df4k)
psp6kres <- data.frame(modified_psp_df6k)
psp8kres <- data.frame(modified_psp_df8k)

colnames(psp1kres) <- c('wm', 'b', 't1k')
colnames(psp2kres) <- c('wm', 'b', 't2k')
colnames(psp4kres) <- c('wm', 'b', 't4k')
colnames(psp6kres) <- c('wm', 'b', 't6k')
colnames(psp8kres) <- c('wm', 'b', 't8k')

des_res <- inner_join(psp1kres, psp2kres, by = c("wm", "b")) 
des_res <- inner_join(des_res, psp4kres, by = c("wm", "b")) 
des_res <- inner_join(des_res, psp6kres, by = c("wm", "b")) 
des_res <- inner_join(des_res, psp8kres, by = c("wm", "b")) 
write.csv(des_res,file='../../../Desktop/ECBS_data/combined_res.csv')
