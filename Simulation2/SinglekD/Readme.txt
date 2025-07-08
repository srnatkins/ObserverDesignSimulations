This folder contains files relating to Simulation 2 but with only one Lipschitz Function studied. Parameters for h, tau, and tau1 are different from what was used in Section 7.2 (moreover T_* is smaller too so that the code could run faster). 
Additionally, only one \Delta(\xi,t) function is being studied in these sets of files. 

Parameters: 
A=[0 1; 0 0];
C = [1 0];
h=3.8;
tau1=0;
tau=0.125;
gap=5;   % this parameter is used to construct a time interval that will have T_*>h+tau+tau1
eps =1;  %epsilon

Some function handles constructed here 
difu.m      Constructs system of differential equations relating to Eqn (30) and (31) from Text
            %Code writes system of 4 differential equations
            % u = [\hat\xi_1,\hat\xi_2,x_1,x_2]^{\top} 
            % RHS of the dynamical system came from (30) and (31) when \delta_a=0 and 
            % \phi(\cdot)=0. 
            % The nonlinearity function used for the simulation is 
            % Del(z)=[kd*sin(z_1(t));0] where z(t) is valued in \R^2
            %term del = [Del(\hat\xi(t)); Del(\hat\xi(t) + x(t))-Del(z(t))] \in R^4
            % note the last two rows of del correspond to \Delta_d(x(t)) as described in (32)
            %A = [0 1; 0 0] will be converted into a 4 by 4 block matrix bigA
            %bigA = A 0 
            %       0 A
            %Within this function rho1=[0;1];
            Inputs: t, u , eps, kd
                    t and u -> symbolic variables needed when calling ode45, 
                    eps     -> epsilon (as described in text)
                    kd      -> Lipschitz constant which will be set based upon function handle Assumption5_1.m
            Output: Differential equation from (30) and (31)

usolve.m   This function solves ode system (30)-(31) as written in difu
            %INPUTS
            % tint  -> time interval tint=[0,T_*] 
            % eps   -> uncertainty parameter
            % kd    -> the lipschitz constant used in \Delta
            % u0    ->  R^2 vector that is the initial condition for (30) and (31) 
                        In main1.m we set u0=[1;0;1;0]; 

            %OUTPUTS
             t     -> vector whose entries correspond to the mesh points used when partitioning tint solving ode system. 
             u     -> a length(t) by 4 matrix with the (i,j) being  u(i,j)=u_j(t_i)
             xihat -> a lenght(t) by 2 matrix representing the first two columns of u
             x     -> length(t) by 2 matrix representing the last two columns of u

ysolve.m        Function handle computes y based upon the spline interpolation of the state 
                Inputs:   tint, eps, u0, kd (same description as above), and tau (time delay) 
                Output:   y (computed based upon equation (29) and then interpolated through spline interpolation)

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

getTstar.m      Function is used to ensure that tint=[0,T_*] is large enough to satisfy Assumption 2.2
                Inputs: tau,tau1,h, gap
                        gap is any positive number depending on how long you want T_* to be                 Output: Tstar, tint, L
                        L = h+tau1+tau
                        Tstar=2*L+gap
                        tint = [0,tstar]

                        

findC0.m        Function handle computes C0=C*[inverse of exp(-A(tau))]. 
                Inputs: C and tau
                        C = [1,0,0,0,0]  
                        tau -> time delay
                Output: C0

findE.m      Function handles computes E and the inverse of E (which is notated as E1). 
              Expression for E is given in Eqn. (8) of text. We use MATLAB's integral function for approximating the integral described in (8)
              Input: h
              Output: E and E1 (inverse of E)

Assumption5_1.m             
            Function returns betastar,Lstar, betabarnew, and kd for simulation 2
            %INPUTS: 
            %       tint        time interval [0,T_*]
            %       tau         delay
            %       h    
            %       tau1
            %%NOTE it is imporitant that inputs are chosen to where
            % T_^*>h+tau1+tau hold as this is one condition needed for Assumption 2.2 
            %OUTPUTS:
            %       betastar    constructed based upon assumption 2
            %                   where beta3 is constructed via Remark 12  
            %       Lstar       the left inverse of betastar           
            %       betabarnew  constructed via equation (67)
            %       kd          set kd=1/(2betabarnew) (must check that kd>0)

            %Based upon sim. 2 computations, betastar is a constant function so
            %betastar is returned as number. Additionally, Lstar=1/betastar is a number 
            %which exists so long as betastar is nonzero.  This is along with 
            % T_*>h+tau1 will make it to where Assumption 2.2 to be satisfied

            %kd>0 chosen such that kd*betabarnew<1 which is a condition for Assumption 5.1. 

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

getmustar.m  Function handle constructs \mu_* (written as mustar in codes). Computation is split into two cases where case 1 relates to Remark 3.7 and Case 2 is otherwise. 
             The later case has mustar computed based upon Eqn (10)

getDs.m      Function computes D_#,D_*, and D_{**} (it returns D_# and D_{**}) based upon eqn. 37
              Inputs: tint,eps,u0,kd,tau,tau1,h
              Outputs: Dsharp, Ddoublestar

getxhatandepshat.m



            
