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
            % t     -> vector whose entries correspond to the mesh points used when partitioning tint
            % solving ode system. 
            % u     -> a length(t) by 4 matrix with the (i,j) being  u(i,j)=u_j(t_i)
            % xihat -> a lenght(t) by 2 matrix representing the first two columns of u
            % x     -> length(t) by 2 matrix representing the last two columns of u



            
