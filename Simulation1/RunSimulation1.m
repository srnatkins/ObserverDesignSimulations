function main
clear all
clf
tint=[0,6]
k = 1
disp('k=1 implies we are using basis functions that are Guassian form')
tau = 0.2
h = 1
tau1 = 1.5
tau2 = 1.7
tau3 = 2.25
eps = [1,2,3]
%x0 =[10;4;0.5;0.25;0]
%x0 = [1,1,1,1,1]';
%x0 = [50;25;1;1;1]
%x0 =[1;20;20;5;10]
x0 = [1;-2;3;-4;5];
%computing solution to x and saving values for plotting purposes
[tplot,x]= xsolve(tint,eps,k,x0);
xplot = transpose(x);
x1plot = xplot(1,:);
x2plot = xplot(2,:);
x3plot = xplot(3,:);
x4plot = xplot(4,:);
x5plot = xplot(5,:);



%% getting beta_** and printing beta_**
[betastar,t,beta11,beta12,beta13] = getbetastar(tint,k,tau,h,tau1,tau2,tau3);

disp('betastar')
disp(betastar)
%computing determinant of betastar
disp('det(betastar)')
det(betastar)
%% computing and printing L_**
L = inv(betastar)
disp('checking if L_** is an inverse of beta_** by multiplying L_** and beta_**')
L*betastar
%% finding mu_**
A = [0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
C = [1 0 0 0 0];

y = ysolve(tint,eps,k,x0,tau);

fun = @(s) expm(transpose(A)*s)*transpose(C)*ppval(y,s);
%% computing inner integral given in (10) at tau1,tau2,and tau3
Int_tau1 = integral(fun,tau1-h,tau1,'ArrayValued',true);
Int_tau2 = integral(fun,tau2-h,tau2,'ArrayValued',true);
Int_tau3 = integral(fun,tau3-h,tau3,'ArrayValued',true);
[E,E1] = findE(h);
%% computing mu(tau1) mu(tau2) and mu(tau3)
mu_tau1 = ppval(y,tau1)-C*E1*expm(-transpose(A)*tau1)*Int_tau1;
mu_tau2 = ppval(y,tau2)-C*E1*expm(-transpose(A)*tau2)*Int_tau2;
mu_tau3 = ppval(y,tau3)-C*E1*expm(-transpose(A)*tau3)*Int_tau3;
disp('mu_** is')
mu_star2 = [mu_tau1,mu_tau2,mu_tau3]

da = @(t) 3*[1*sin(t); 2*sin(2*t); 3*sin(3*t); 4*sin(4*t);5*sin(5*t)];
db = @(t) 3*sin(t);
[tdelta, dc,dd,de,df] = getds(tint,da,db,tau,h,tau1,tau2,tau3);
checkeps = mu_star2*transpose(L)+transpose(df)*(transpose(L))

 xu = NaN(5,length(t));
 for i = 1:length(t)
     beta1 = [beta11(:,i),beta12(:,i),beta13(:,i)];
     xu(:,i) = E1*expm(-transpose(A)*t(i))*integral(fun,t(i)-h,t(i),'ArrayValued',true)+beta1*L*transpose(mu_star2)+dd(:,i)+beta1*L*df;
 end
xu1 = xu(1,:);
xu2 = xu(2,:);
xu3 = xu(3,:);
xu4 = xu(4,:);
xu5 = xu(5,:);

figure(1) 
plot(tplot,x1plot,'b-', t, xu1, 'b-.','LineWidth', 2);
hold on 
plot(tplot,x2plot,'r-', t, xu2, 'r-.','LineWidth', 2);
hold on 
plot(tplot,x3plot,'g-',t, xu3, 'g-.','LineWidth', 2);
hold on 
plot(tplot,x4plot,'c-',t, xu4, 'c-.','LineWidth',2);
hold on 
plot(tplot,x5plot,'k-',t,xu5,'k-.','LineWidth',2);
xline(h+tau3+tau)
legend('$x_1(t)$','$x_{u1}(t)$','$x_2(t)$','$x_{u2}(t)$','$x_3(t)$','$x_{u3}(t)$','$x_4(t)$','$x_{u4}(t)$','$x_5(t)$','$x_{u5}(t)$', '$h+\max{\tau_i}+\tau$','Interpreter', 'latex', 'FontSize', 11)
ylim([-115,125]) 
ax = gca;
ax.FontSize = 20;
yplot = zeros(size(tplot));
for i = 1: length(tplot)
    yplot(i)=ppval(y,tplot(i));
end
figure(2)
plot(tplot, x1plot, 'r-', tplot,yplot,'b-.','LineWidth',2);
xline(h+tau3+tau);
legend('$x_1(t)$', '$y(t)$', '$h+\max{\tau_i}+\tau$','Interpreter','latex','FontSize',18)
ax=gca;
ax.FontSize=20;

figure(3)
plot(tplot,x1plot,'r-', tplot,yplot,'b-.','LineWidth',2);
xline(tau);
legend('$x_1(t)$','$y(t)$', '$\tau$','Interpreter','latex','FontSize',18)
ax=gca;
ax.FontSize=20;


end