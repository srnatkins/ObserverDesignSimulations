function [tvec,xhat,epshat] = getxhatandepshat(tint,eps,u0,kd,Lstar, tau,tau1,h)
[E,E1] = findE(h);
A=[0 1; 0 0];
C = [1 0];
%call beta1solve to get beta1 
[tvec,beta1] =beta1solve(tint,tau,h);

%call ysolve to get y(s)
y = ysolve(tint,eps,u0,kd,tau);
% call getDs to get Dsharp and Ddoublestar
[Dsharp,Ddoublestar]= getDs(tint,eps,u0,kd,tau,tau1,h);
% call getmustar
mustar = getmustar(tint,eps,u0,tau,h,tau1,kd); %returns a vector whose entries pertain to mu*(tvec(i))
func1 = @(s) expm(transpose(A)*s)*transpose(C)*ppval(y,s);
term1 = zeros(2,length(tvec));
term2 =zeros(2,length(tvec));


for i = 1: length(tvec)
    term2(:,i) = Dsharp(tvec(i))+beta1*Lstar*(mustar(i)+Ddoublestar(tvec(i)));
    term1(:,i) = E1*expm(-transpose(A)*tvec(i))*integral(func1,tvec(i)-h,tvec(i), "ArrayValued",true);
    Dstarstar(i) = Ddoublestar(tvec(i));
end
xhat = term1+term2;
epsstar = mustar+Dstarstar;
epshat = epsstar*Lstar;

end