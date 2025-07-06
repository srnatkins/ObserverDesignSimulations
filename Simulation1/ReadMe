This repository gives the files needed to compute the simulations performed in Section 7.1 of the article titled "Fixed and almost fixed time estimation under unknown parameters in measurements" (Authors: Summer Atkins, Michael Malisoff, and Frederic Mazenc). 
The simulations provided are all MATLAB files. The simulations relate to a specific dynamical system of a single-link robotic manipulator coupled to a DC motor, which is of the form Eqn. (3) from the text.
Be advised that these codes are written specifically for this problem.
Within these files we take advantage of $\gamma_{ai}$'s being set to zero (this has an effect on how the $\gammma_i$'s is computed). 
Additionally, h, $\tau$, $\tau_1$, $\tau_2$, and $\tau_3$ were chosen to where Remark 3.7 can be leveraged.
Note: Some objects (variables, functions, parameters) used within the article were renoteded with in these codes for shorthand purposes and also because some symbols like "*" cannot be part of the name of a variable in MATLAB.
Some examples include the following: 
$\delta_a$, $\delta_b$,...,$\delta_f$ were renotated as "da", "db", ..., "df" respectively, 
Vector $\epsilon$ was renotated as "eps",  
$\lambda_i$ was renotated as the ith entry of vector "lam",
$\beta_{**}$ (see Remark 3.7) is renotated as "betastar",
$\mu_{**}$ (see Remark 3.7) is renotated as "mustar" 
$\L_{∗∗}$ (see Remark 3.7) is renotated as "L" 
$E^{-1}$ (inverse of E) is renotated as E1
$\rho_i(t)$ is renotated as rhoi whose 

Addtionally: dimensional parameters p, n, and q are fixed and not explicitly defined in the files, 
but rather implicitly defined through settings for certain objects within the file. For this example, p=1, n=5, and q=3. 
Also some objects within the article (ie matrix A, C, $\lambda_i$, w, $\rho_i$) maybe defined in multiple places within the files. 
  

There are many m files within this repository.
% RunSimulation1.m  -> Main file used for running simulations

  difx.m       Function handle that creates the dynamical system for the Motor example to be solved via ode45 
                Inputs: t,x, eps, k  
                              t and x are symbolic variables representing independent variable t and state variable x
                              eps  is a vector of dimension q (for this example q is 3) of unknown real values
                              k    is a case based parameter where users can choose between the following basis function forms in constructing function $\rho_i$
                                   k = 1  -> Guassian form                 f(t) =exp(-t)  this input is written in RunSimulation1.m and set as 1.
                                   k = 2  -> Multiquadratic form           f(t) =sqrt(1+t)
                                   k = 3  -> Inverse Multiquadratic form   f(t) = 1/sqrt(1+t)
                                   k = 4  -> Cauchy form                   f(t) =  1/(1+t)        
              Output: The differential equation as described in (3) but with settings based upon the motorexample as described in section 7.1 that is
                        C=[1 0 0 0 0], 
                        A as described in (61) but with all parameters used within those matrix entries set to 1, 
                        lam=[1 2 3], w=[1,2,3],   (note w(i) and lam(i) is matlab notation for the ith entry of vectors w and lam)
                        Basis function used for defining $\rho_i$ is Guassian form f(t) = exp(-t)  (that is k=1 in our simulations)
                        rho1 is stored as a symbolic vectored valued function valued in R^n (n=5 in this case) in which each component is set as 
                              f( (t-lam(1))^2/(2w(1)^2) )  where f(t) = exp(-t)
                        Similarly rho2 and rho3 are symbolic function valued in R^n (n=5) with each component is 
                        set as f( (t-lam(2))^2/(2w(2)^2) ) and f( (t-lam(3))^2/(2w(3)^2) ) respectively. 
                        $\delta_a$ which is notated as da and stored as an   n vector valued function whose components are defined as 
                            The ith entry of da = 3*i*sin(i*t)  for i = 1,...,n (n=5)   
                          
  xsolve.m     Function handle that uses ode45 to solve the differential equation generated from difx.m
          Notes: Line 4 are ode45 setting adjustment relating to relative tolerance, absolute tolerance, and maximum stepsize. 
                If parameter changes to the files lead to an error in running xsolve, one may consider adjusting the tolerances.
                Inputs: t, eps, k, x0
                      t=[t0,tf]  -> a 2 vector representing the intial time t0 and final time tf for which the problem is being studied
                      eps        -> q vector of unknown real values used in the dynamical system
                      k          -> case based parameter (see above for description)
                      x0         -> initial condition used for state system
                Output: [t,x] 
                      t          -> vector whose entries represent the mesh points used for partitioning interval [t0,tf]
                      x          -> Array that generates the approximated value of the states at the mesh points

                  
              
interpolate_x.m Function handle takes numerical sol. from xsolve.m and performs spline interpolation 
                Inputs:   t, eps, k, x0 (same description as above)
                Output:   x interpolated through spline interpoltation

ysolve.m        Function handle computes y based upon the spline interpolation of the state 
                Inputs:   t, eps, k, x0 (same description as above), and tau (time delay) 
                Output:   y (computed based upon equation (3) and then interpolated through spline interpolation)

findC0.m        Function handle computes C0=C*[inverse of exp(-A(tau))]. 
                Inputs: C and tau
                        C = [1,0,0,0,0]  
                        tau -> time delay
                Output: C0

findE.m      Function handles computes E and the inverse of E (which is notated as E1). 
              Expression for E is given in Eqn. (8) of text. We use MATLAB's integral function for approximating the integral described in (8)
              Input: h
              Output: E and E1 (inverse of E)

gamsolve.m    Function handle computes $\gammma_{i}$ for i=1,...,q (q=3 in this example) by using Eqn (5) from text. 
              MATLAB's integral function (with 'ArrayValued' settings turned on) is used for approximating integrals given here. 
              $\gamma_{bi}$ and then leverages the fact that $\gammma_{ai}=0$ in the example to return that $\gamma_{i}=\gammma_{bi}$. 
              Naturally this code can be modified to where nonzero values of $\gamma_{ai}$ can be stored and thus return $\gamma_{i}=\gamma_{ai}+\gamma_{bi}$ 
              Inputs: k and tau (as described in prevous function handles but some more clarification is given about k below)
                      k is the case based parameter for determining what type of basis function is used for constructing $\rho_{i}$
              Outputs: gam1 gam2 and gam3  as described in Eqn (5)

beta1solve.m    Function handle used for approximated \beta_{1i}(t) for i = 1,...,q (q=3) as defined in the first equation of Eqn. (8)
                MATLAB's integral function (with 'ArrayValued' settings turned on) is used in approximating intergrals used in Eqn. (8)
                Inputs:  tint, k, tau, h
                    tint=[t0,tf]   -> a 2 vector representing the intial time t0 and final time tf for which the problem is being studied  
                    k              -> cased based parameter (see detailed description above)
                    tau            -> time delay
                    h              -> constant (should be the case that h>tau) used in Eqn. 8 of text
                 Outputs: t, beta11, beta12, beta13
                          t vector of length N whose entries are the mesh points used in partitioning time interval  [t0, tf]. 
                            NOTE: The partition being used for the beta functions are equally spaced of mesh size ds (currently set to being 0.05). 
                                  The purpose of setting ds=0.05 is to ensure choices of constants tau1, tau2, and tau3 are part of the partition (more on this in getbetastar.m) 
                          beta11, beta12 and beta13 are based on the first equation in (8). 
                          Approximation of the functions values of beta11,beta12,beta13 at the mesh time points are stored. 

beta3solve.m    Function handle used for approximating \beta_{3i}(t) for i = 1,..,q based upon the second to last equation in (8). 
                \beta_{2i} is used within this code but not needed elsewhere so \beta_{21} is not a saved variable. 
                Inputs: tint, k, tau , h tint=[t0,tf] and k, tau, and h are as defined before. 
                Outputs: time vector t, Vectors Beta31, Beta32, and Beta33, and matrices beta11, beta12, beta13
                        t                     -> partitioning of tint
                        Beta31,Beta32,Beta33  -> vectors of the same size as t but with entries representing approximations of fucntion $\beta_{31}$ as the mesh time points. 
                        beta11,beta12,beta13  -> see description above.

interpolatedbeta1.m   Function handle interpolated beta1i through MATLAB's spline function 
                      Inputs: t, beta11, beta12, and beta13 
                      Outputs Interpolating functions Beta11 Beta12 Beta13 through spline interpolation
                      Currently this function is NOT being used. 

getbetastar.m         Computing $\beta_{**}$. That is the code uses Remark 3.7 instead of using eqn (9) to approximate $\beta_*$
                      Inputs: tint, k, tau, h, tau1, tau2, tau3
                              tau, h, tau1, tau2, and tau3 are the constants that are picked to where there is some positive integer N (in this case N=3) where
                                    \tau_i \geq h+tau for i = 1,...,N  (here \geq is latex command meaning greater than or equalt to) 
                              This condition is needed in order to use Remark 3.7. Within this code \beta_{31}, \beta_{32}, and \beta_{33} are computed via beta3solve then  
                              satisfy taui is greater than or equal to h+tau for all i=1,2,3. 
                      Outputs: betastar, t, beta11, beta12, beta13  (t, beta11, beta12, beta13 are same as described in beta1solve.m)
                      $\beta_{**}$ is needed for computing $L_{**}$

getds                 Computes the values of $\delta_c$, $\delta_d$, $\delta_e$, and $\delta_f$
                      From the article $\delta_c$ is given in Eqn (7), and $\detal_d$, $\delta_{e}$,and $\delta_f$ are given in Eqn. (11).
                      Inputs: tint, da, db, tau, h, tau1, tau2, tau3
                              da is \delta_a which is set to  3*[1*sin(t); 2*sin(2*t); 3*sin(3*t); 4*sin(4*t);5*sin(5*t)];
                              db is \delta_b which is set to 3*sin(t)      
                      Output: t, dc,dd,de,df

Runsimulation1.m     Generates what x and epsilon should be based upon Theorem 2.3
                     Within this code various inputs needed for the function handles above are written. 
                     Then function handles are called. Then x is approximated based upon Eqn (12) from Theorem 2.3. 
                     Comparison plots between x and the estimator of x are then generated. 



                      

                          
                        
                
                          
                    
                    
      
