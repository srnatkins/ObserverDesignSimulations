function [t,beta11,beta12,beta13]= beta1solve(tint,k,tau,h)
%defining rhoi(t)
lam =[1,2,3];
w =[1,2,3];
if k == 1 %Guassian form
   f = @(t) exp(-t);
elseif k == 2 %multiquadratic
    f = @(t) sqrt(1+t);
elseif k == 3 %inverse multiquadratic
    f = @(t) (1+t)^(-1/2);
elseif k == 4 %cauchy form
    f = @(t) (1+t)^(-1);
end
rho1 = @(t) [f( (t-lam(1))^2/( 2*w(1)^2 ) );f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) )];
rho2 = @(t) [f( (t-lam(2))^2/( 2*w(2)^2 ) );f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) )];
rho3 = @(t) [f( (t-lam(3))^2/( 2*w(3)^2 ) );f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) )];

%defining A and C
A = [0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
C = [1 0 0 0 0];
% creating time interval
dt=0.05;
t = tint(1):dt:tint(end);
%intializing functions beta1i---i corresponds to rhoi(t)
beta11 = zeros(5,length(t));
beta12 = zeros(5,length(t));
beta13 = zeros(5,length(t));
%finding E and Einver
[E,E1] = findE(h);
%definining beta1i(t)
for i =1:length(t)
    fun1 = @(m) expm(-A*m)*rho1(m);
    fun2 = @(m) expm(-A*m)*rho2(m);
    fun3 = @(m) expm(-A*m)*rho3(m);
    innerint1 = @(s) integral(fun1,s-tau,t(i),'ArrayValued',true);
    innerint2 = @(s) integral(fun2,s-tau,t(i),'ArrayValued',true);
    innerint3 = @(s) integral(fun3,s-tau,t(i),'ArrayValued',true);
    Int1 = @(ell) expm(transpose(A)*ell)*transpose(C)*C*expm(A*ell)*innerint1(ell);
    Int2 = @(ell) expm(transpose(A)*ell)*transpose(C)*C*expm(A*ell)*innerint2(ell);
    Int3 = @(ell) expm(transpose(A)*ell)*transpose(C)*C*expm(A*ell)*innerint3(ell);
    beta11(:,i) = E1*expm(-transpose(A)*t(i))*integral(Int1,t(i)-h,t(i),'ArrayValued',true);
    beta12(:,i) = E1*expm(-transpose(A)*t(i))*integral(Int2,t(i)-h,t(i),'ArrayValued',true);
    beta13(:,i) = E1*expm(-transpose(A)*t(i))*integral(Int3,t(i)-h,t(i),'ArrayValued',true);
end
%old code
% for i = 1:length(t)
%     %create a new time vector that partitions [t(i)-h,t(i)]
%     ts = t(i)-h:1.e-3:t(i);
%     %initialize inner integral see L_5i in equation (58)
%     innerint1 = zeros(5,length(ts));
%     innerint2 = zeros(5,length(ts));
%     innerint3 = zeros(5,length(ts));
%     for j = 1: length(ts)
%         innerfun1 = @(m) expm(-A*m)*rho1(m);
%         innerfun2 = @(m) expm(-A*m)*rho2(m);
%         innerfun3 = @(m) expm(-A*m)*rho3(m);
%         innerint1(:,j) = integral(innerfun1,ts(j)-tau,t(i),'ArrayValued',true);
%         innerint2(:,j) = integral(innerfun2,ts(j)-tau,t(i),'ArrayValued',true);
%         innerint3(:,j) = integral(innerfun3,ts(j)-tau,t(i),'ArrayValued',true);
%     end
%     innerInt1 = spline(ts,innerint1)
%     innerInt2 = spline(ts,innerint2);
%     innerInt3 = spline(ts,innerint3);
%     funb1 = @(s) expm(transpose(A)*(s-t(i)))*transpose(C)*C*expm(A*s)*ppval(innerInt1,s);
%     funb2 = @(s) expm(transpose(A)*(s-t(i)))*transpose(C)*C*expm(A*s)*ppval(innerInt2,s);
%     funb3 = @(s) expm(transpose(A)*(s-t(i)))*transpose(C)*C*expm(A*s)*ppval(innerInt3,s);
%     beta11(:,i) = E1*integral(funb1,t(i)-h,t(i),'ArrayValued',true);
%     beta12(:,i) = E1*integral(funb2,t(i)-h,t(i),'ArrayValued',true);
%     beta13(:,i) = E1*integral(funb3,t(i)-h,t(i),'ArrayValued',true);
% end

end