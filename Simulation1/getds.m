function [t, dc,dd,de,df] = getds(tint,da,db,tau,h,tau1,tau2,tau3)
A = [0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
C = [1 0 0 0 0];
Csh = transpose(C)*C;
[E,E1] = findE(h);
fun1 = @(m) expm(-A*m)*da(m);
dc =@(t) db(t)-C*expm(A*t)*integral(fun1,t-tau,t,'ArrayValued',true);
%computing delta_d
dt = 0.05;
t = tint(1):dt:tint(end);
istar1 = round(tau1/dt)+1;
istar2 = round(tau2/dt)+1;
istar3 = round(tau3/dt)+1;
dd = NaN(5,length(t));
de = NaN(1,length(t));



for i = 1:length(t)
fun2 = @(s) expm(transpose(A)*s)*transpose(C)*db(s);
ddterm1 = -E1*expm(-transpose(A)*t(i))*integral(fun2,t(i)-h,t(i),'ArrayValued',true);
fun3 = @(s) integral(fun1,s-tau,t(i),'ArrayValued',true);
fun4 = @(s) expm(transpose(A)*s)*transpose(C)*C*expm(A*s)*fun3(s);
ddterm2 = E1*expm(-transpose(A)*t(i))*integral(fun4,t(i)-h,t(i),'ArrayValued',true);
dd(:,i) = ddterm1+ddterm2;
 de(i) = dc(t(i))+C*dd(:,i);
end
df = -[de(istar1);de(istar2);de(istar3)];

end