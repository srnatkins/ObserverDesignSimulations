function [t,x] =xsolve(tint,eps,x0)
t0 = tint(1);
tf = tint(end);
options = odeset('RelTol',1e-8,'AbsTol',1e-9, 'MaxStep',1.e-2);
%options = odeset('RelTol',1e-8,'AbsTol',1e-9);
[t,x] =ode45(@(t,x)difx(t,x,eps),[t0,tf],x0,options);
end
