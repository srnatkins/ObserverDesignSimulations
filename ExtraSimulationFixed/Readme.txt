These files provide illustrations on how variation of parameters based methods can be applied to problem (3) but with \delta_a and \delta_b and all of the \gamma_{ai}'s set to 0 (See Remark 6.3). 
The example given relates to the one describe in Section 7.2 but without any \Delta function. 

Inputs in main.m
tint=[0,10];
tau = 0.125;
h = 3.8;
tau1 = 0;
A=[0 1; 0 0];
C=[1,0];
x0 = [1;1];
eps = 2;

difx.m       Function handle that creates the dynamical system for extrasimulation 
                Inputs: t,x, eps  
                              t and x are symbolic variables representing independent variable t and state variable x
                              eps  is a vector of dimension q (for this example q is 1) of unknown real values
                Output: dx = Ax+eps*rho1   % where rho1=[0;1];

 xsolve.m     Function handle that uses ode45 to solve the differential equation generated from difx.m
                Inputs: t, eps, x0
                      t=[t0,tf]  -> a 2 vector representing the intial time t0 and final time tf for which the problem is being studied
                      eps        -> q vector of unknown real values used in the dynamical system
                      x0         -> initial condition used for state system
                Output: [t,x] 
                      t          -> vector whose entries represent the mesh points used for partitioning interval [t0,tf]
                      x          -> Array that generates the approximated value of the states at the mesh points

interpolate_x.m Function handle takes numerical sol. from xsolve.m and performs spline interpolation 
                Inputs:   t, eps, k, x0 (same description as above)
                Output:   x interpolated through spline interpoltation

ysolve.m        Function handle computes y based upon the spline interpolation of the state 
                Inputs:   t, eps, x0 (same description as above), and tau (time delay) 
                Output:   y (computed based upon equation (3) and then interpolated through spline interpolation)

difLi          Function handle generates the delay differential equation presented in Remark 6.3 eqn: (59)
                Inputs: t, y, yL, tau
                        t,y,yL are all symbolic variables yL is used to represent the lag term that will have delay and tau is a constant representing the delay
                Output: the delayed differential equation

yhist.m        stores the history function used for the delay differential equation described in difLi (yhist =[0,1,0,1]^T)

Lsolve.m      Function handle used for numerically solving dynamical system generated from difLi.  
                 Input: vector t=[t0,tf] and time delay tau
              When calling dde23 settings for the relative tolerance, absolute tolerance, and maximum stepsize are adjusted via ddeset command.
              Additionally, @yhist is an input for dde23 
                  Output: t1, L11, L21, sol1
                          t1    -> time vector whose entries correspond to the nodes used in partitioning time interval [t0,tf]
                          L11   -> Approximation of the solution to the  first dde equation given in (59) at the time values given in time entries to t1
                          L21   -> Approximation of the solution to the  second dde equation given in (59) at the time values given in time entries to t1
                          sol1  -> L11 and L21 combined (makes things easier for the function handled described below)
             
interpolatL.m Function handle calls Lsolve and then interpolates the numerical approximation of the delayed differential equation via MATLAB's spline function. 
                Input: t=[t0,tf] and time delay tau
                Output: interpolated functions L11 and L21 


findC0.m      Function handle computes C0=C*[inverse of exp(-A(tau))]. 
                Inputs: C and tau
                        C = [1,0]  
                        tau -> time delay
                Output: C0

findE.m      Function handles computes E and the inverse of E (which is notated as E1). 
              Expression for E is given in Eqn. (8) of text. We use MATLAB's integral function for approximating the integral described in (8)
              Input: h
              Output: E and E1 (inverse of E)

getbeta.m    Function handle used for generated beta_{*} 
             First interpolateL is called. Then computations for \beta_{*} (notated as betastar within the code) is split into 2 cases. 
             Case 1: The Inputs h, tau, and tau1 are chosen to where tau1 is greater than or equal to tau+h (ie the case in which Remark 3.7 is being used). 
                     For this case \beta_{*}=\beta_{**}=\beta_3(tau1) where the expression of \beta_3(t) is given in Eqn (57). 
                     Note: this case should also ensure that \beta_{**} has a left inverse but this is checked later within the text. 
             Case 2: Otherwise \beta_{31} computed based upon equation (57), and then betastar is set based upon Eqn (9) of text that is betastar=beta_{31}(t-tau1).
                     Note in Case 2, betastar requires timer interval [t0,tf] to be partitioned. In this case time is partitioned by subdividing [t0,tf] into a fixed number of equally spaced nodes with mesh                      size being 0.1. Consequently, betastar is stored as an array with components relating to approximations of \beta_{*} valued at the nodes.  
             Inputs: t=[t0,tf], tau, tau1, and h
             Output: betastar

beta1solve.m Function handle used for computing beta1 from Eqn 8  
             Inputs: tint=[t0,tf], tau, h
             Output: time vector t and beta1

getmustar.m  Function handle constructs \mu_* (written as mustar in codes). Computation is splint into two cases where case 1 relates to Remark 3.7 and Case 2 is otherwise. 
             The later case has mustar computed based upon Eqn (10)

getMS.m      Function isn't used but would related to computations for what would be needed in Assumption 5.1/Theorem 5.2. This is more so important with regard to Problems of the form Eqn. (29)

main.m       Script file that has the input settings as described at the top of the Readme file and calls the function handles above to generate an estimator x_u and \epsilon_u via Eqn (13) (note this is the same as equation 12 since \delta_a,\delta_b, gamma_{a1} are set to 0). 
             File prints values of \epsilon_u (which is labled as epscheck) and plots xu which is compared to the true solution x to demonstrate fixed time estimation. 
             

                            


                
      

