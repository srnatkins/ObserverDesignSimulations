function main
clear
clear fig
x0 = [1;1];
eps = 2;
tint =[0,5];
[t,x] =xsolve(tint,eps,x0);
tsave=t;
xplot1= x(:,1);
xplot2 = x(:,2);

A=[0 1; 0 0];
C=[1,0];


tau = 0.25;
h=1.2;

tau1 = 0;

betastar = getbeta(t,tau,tau1,h);
Lstar= 1./betastar;



[t,beta1] =beta1solve(tint,tau,h);
% format long
% beta1
% format short

mustar = getmustar(tint,eps,x0,tau,h,tau1);


epscheck= mustar.*Lstar;
epscheck=epscheck(end)
[E,E1] = findE(h);
 y = ysolve(tint,eps,x0,tau);
xu = NaN(2,length(t));
t = tint(1): 0.1:tint(end);
fun = @(s) expm(transpose(A)*s)*transpose(C)*ppval(y,s);
for i = 1:length(t)
    xu(:,i) = E1*expm(-transpose(A)*t(i))*integral(fun,t(i)-h,t(i),'ArrayValued',true)+beta1(:,i)*Lstar(i)*mustar(i);
end
xu1 =xu(1,:);
xu2 =xu(2,:);
figure(1)
plot(tsave,xplot1, 'b-', t, xu1,'b-.', 'LineWidth',2);
hold on
plot(tsave,xplot2,'r-', t, xu2,'r-.', 'LineWidth', 2);
hold on
xline(h+tau1+tau)
legend('$x_1(t)$','$x_{u1}(t)$','$x_2(t)$','$x_{u2}(t)$','$h+\max{\tau_i}+\tau$','Interpreter', 'latex', 'FontSize', 12)
M1=getMs(tint,tau,h,tau1)

end