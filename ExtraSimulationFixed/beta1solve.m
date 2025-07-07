function [t,beta1] =beta1solve(tint,tau,h)
rho1 = [0;1];
A = [0,1;0,0];
C = [1,0];
dt = 0.1;
t = tint(1):dt:tint(end);
beta1 = NaN(2,length(t));
[E,E1] = findE(h);
for i = 1:length(t) 
    fun1 = @(m) expm(-A*m)*rho1;
    innerint1 = @(s) integral(fun1,s-tau,t(i),'ArrayValued',true);
    Int1 = @(ell) expm(transpose(A)*ell)*transpose(C)*C*expm(A*ell)*innerint1(ell);
    beta1(:,i) = E1*expm(-transpose(A)*t(i))*integral(Int1,t(i)-h,t(i),'ArrayValued',true);
end
end