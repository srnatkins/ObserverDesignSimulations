function M1=getMs(tint,tau,h,tau1)
A =[0 1; 0 0];
C = [1 0];
Csharp= transpose(C)*C;
[E,E1] = findE(h);  %E1 is the inverse of E
[t,beta1] =beta1solve(tint,tau,h);
% disp(beta1)
betastar = getbeta(tint,tau,tau1,h);
Lstar = 1./betastar;
M1 = beta1.*Lstar; %will be constant for t(i)>tau+h+tau1
M1=M1(:,end);
fun1 = @(l) norm(M1*C*expm(-A*l));
term1 = integral(fun1,-tau,0,'ArrayValued',true);

fun2 = @(s,m) norm(M1*C*E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
fun3 = @(s,m) norm(E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
term2 = integral(@(s) integral(@(m)fun2(s,m), s-tau,0,'ArrayValued',true), -h,0,'ArrayValued',true)
term3 = integral(@(s) integral(@(m)fun3(s,m), s-tau,0,'ArrayValued',true),-h,0,'ArrayValued',true)

fun4 = @(s,m) norm((eye(2)-M1*C)*E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
term4 = integral(@(s) integral(@(m)fun4(s,m), s-tau,0,'ArrayValued',true),-h,0,'ArrayValued',true)

betabar1 = term1+term2+term3
betabar2=term1+term4
check=term1+term2-term3


end