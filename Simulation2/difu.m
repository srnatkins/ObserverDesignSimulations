%Code writes system of 4 differential equations
% u = [\hat\xi_1,\hat\xi_2,x_1,x_2]^{\top} 
% RHS of the dynamical system came from (30) and (31) when \delta_a=0 and 
% \phi(\cdot)=0. 
% The nonlinearity function used for the simulation is 
% Del(z)=[kd*sin(z_1(t));0] where z(t) is valued in \R^2
%term del = [Del(\hat\xi(t)); Del(\hat\xi(t) + x(t))-Del(z(t))] \in R^4
% note the last to rows of del correspond to \Delta_d(x(t)) as described in
% (32)
%A = [0 1; 0 0] will be converted into a 4 by 4 block matrix bigA
%bigA  A 0 
%      0 A
%preset inputs rho1=[0;1];
function du = difu(t,u,eps,kd)
    rho1 = [0;1];
    A = [0,1;0,0];
    bigA = [A,zeros(2);zeros(2),A];
    del = [kd*sin(u(1));0;kd*sin(u(3)+u(1))-kd*sin(u(1));0];
    du = bigA*u+eps*[0;0;rho1]+del;
end