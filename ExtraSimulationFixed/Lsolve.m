function [t1, L11, L21,sol1] = Lsolve(t,tau)

lags = tau;
t0 = t(1);
tf = t(end);
options = ddeset('RelTol',1e-8,'AbsTol',1e-9,'MaxStep',1.e-2);
sol1 = dde23(@(t,y,yL)difLi(t,y,yL,tau),lags,@yhist,[t0,tf],options);
t1 = sol1.x;
Lnew = sol1.y;
L11 = Lnew(1:2,:);
L21 = Lnew(3:4,:);



end