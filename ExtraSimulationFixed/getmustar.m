function mustar = getmustar(tint,eps,x0,tau,h,tau1)
y = ysolve(tint,eps,x0,tau);
[E,E1] = findE(h);
C=[1,0];
A = [0 1; 0 0];
fun = @(s) expm(transpose(A)*s)*transpose(C)*ppval(y,s);
if tau1 >h+tau || tau1==h+tau
    Int_tau1 = integral(fun,tau1-h,tau1,'ArrayValued',true);
    mustar = ppval(y,tau1)-C*E1*expm(-transpose(A)*tau1)*Int_tau1;
else
    t = 0:0.1:tint(end);
    mustarvec=zeros(size(t));
    for i = 1: length(t)
        Int= integral(fun,t(i)-tau1-h,t(i)-tau1,"ArrayValued",true);
        mustarvec(i)=ppval(y,t(i)-tau1)-C*E1*expm(-transpose(A)*(t(i)-tau1))*Int;
    end
    %mustarvec(i) is constant when t(i)>h+tau1+tau
    mustar =mustarvec;
end
end