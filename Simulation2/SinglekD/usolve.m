%% This function solves ode system (30)-(31) as written in difu
%INPUTS
% tint  -> time interval 
% eps   -> uncertainty parameter
% kd    -> the lipschitz constant used in \Delta
% u0    ->  R^2 vector that is the initial condition for (30) and (31) 

%OUTPUTS
% t     -> vector whose entries correspond to the mesh points used when
% solving ode system. 
% u     -> length(t) by 4 matrix whose u(i,j)=u_j(t_i)
% xihat -> lenght(t) by 2 matrix representing the first two columns of u
% x     -> length(t) by 2 matrix representing the last two columns of u
%
function [t,u,xihat,x] = usolve(tint,eps,u0,kd)
t0 = tint(1);
tf = tint(end);
options = odeset('RelTol',1e-8,'AbsTol',1e-9, 'MaxStep',1.e-2);
%u0=[x0;x0]; %made the initial conditions for xi the same as x
[t,u] = ode45(@(t,u)difu(t,u,eps,kd),[t0,tf],u0,options);
xihat = u(:,1:2);
x = u(:,3:4);
end