function [t,x]= xsolve(t,eps,k,x0)
t0 = t(1);
tf = t(end);
options = odeset('RelTol',1e-8,'AbsTol',1e-9, 'MaxStep',1.e-2);
[t,x] =ode45(@(t,x)difx(t,x,eps,k),[t0,tf],x0,options);
end