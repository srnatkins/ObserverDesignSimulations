function main
A=[0 1; 0 0];
C = [1 0];
h=3.8;
tau1=0;
tau=0.125;
gap=5;
eps =1;
u0=[1;0;1;0];
[E,E1]=findE(h);
[Tstar,tint,L]=getTstar(tau,tau1,h,gap);
Ldouble = 2*L;
fprintf('L=h+tau1+tau=%f.\n',L);
fprintf('T_star is %f.\n',Tstar);
fprintf('2L is %f.\n', Ldouble);
fprintf('T_star>2L since T*-2L=%f.\n',gap);
[betastar, Lstar, betabarnew,kd] = Assumption5_1(tint,tau,h,tau1);

[tvec,xhat,epshat] = getxhatandepshat(tint,eps,u0,kd,Lstar, tau,tau1,h);
xhat1 = xhat(1,:);
xhat2 = xhat(2,:);
 Trueeps=eps*ones(size(tvec));

[t,u,xihat,x] = usolve(tint,eps,u0,kd);
x=transpose(x);
x1plot=x(1,:);
x2plot=x(2,:);

figure(1)
plot(t,x1plot,'b-','LineWidth',1);
hold on 
plot(tvec, xhat1, 'b-.','LineWidth', 2);
hold on 
plot(t, x2plot,'r-','LineWidth',1);
hold on 
plot(tvec,xhat2,'r--','LineWidth',2);
hold on 
xline(Ldouble)
legend('$x_1(t)$', '$\hat{x}_1(t)$', '$x_2(t)$', '$\hat{x}_2(t)$', '2L','Interpreter', 'latex', 'FontSize', 11)
xlim([0,Tstar])
ax = gca;
ax.FontSize=20;

figure(2)
plot(tvec,epshat, 'k-.','LineWidth',2);
hold on
plot(tvec, Trueeps,'g--','LineWidth',2)
hold on
xline(Ldouble)
legend('$\epsilon_{*}(t)^{\top}L_{*}^{\top}$','$\epsilon$','$2L$','Interpreter', 'latex', 'FontSize', 11 )
xlim([0,Tstar])
ax = gca;
ax.FontSize=20;

figure(3)
plot(t,x1plot,'b-','LineWidth',1);
hold on 
plot(tvec, xhat1, 'b-.','LineWidth', 2);
hold on 
xlim([Ldouble,Tstar])
legend('$x_1(t)$', '$\hat{x}_1(t)$', 'Interpreter', 'latex', 'FontSize', 11)
ax = gca;
ax.FontSize=20;

figure(4)
plot(t, x2plot,'r-','LineWidth',1);
hold on 
plot(tvec,xhat2,'r--','LineWidth',2);
hold on 
xlim([Ldouble,Tstar])
legend('$x_2(t)$', '$\hat{x}_2(t)$', 'Interpreter', 'latex', 'FontSize', 11)
ax = gca;
ax.FontSize=20;

figure(5)
plot(t,x1plot,'b-','LineWidth',1);
hold on 
plot(tvec, xhat1, 'b-.','LineWidth', 2);
hold on 
xlim([10,12])
legend('$x_1(t)$', '$\hat{x}_1(t)$', 'Interpreter', 'latex', 'FontSize', 11)
ax = gca;
ax.FontSize=20;

figure(6)
plot(t,x2plot,'r-','LineWidth',1);
hold on 
plot(tvec, xhat2, 'r--','LineWidth', 2);
hold on 
xlim([10,12])
legend('$x_2(t)$', '$\hat{x}_2(t)$', 'Interpreter', 'latex', 'FontSize', 11)
ax = gca;
ax.FontSize=20;

figure(7)
plot(tvec,epshat, 'k-.','LineWidth',2);
hold on
plot(tvec, Trueeps,'g--','LineWidth',2)
hold on
xline(Ldouble)
legend('$\epsilon_{*}(t)^{\top}L_{*}^{\top}$','$\epsilon$','$2L$','Interpreter', 'latex', 'FontSize', 11 )
ax = gca;
ax.FontSize=20;

figure(8)
plot(t,x2plot,'r-','LineWidth',1);
hold on 
plot(tvec, xhat2, 'r--','LineWidth', 2);
hold on 
legend('$x_2(t)$', '$\hat{x}_2(t)$', 'Interpreter', 'latex', 'FontSize', 11)
ax = gca;
ax.FontSize=20;

[xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd);
xtrue = zeros(2,length(tvec));
epsdiff = zeros(size(tvec));
xdiff = zeros(size(tvec));
for i = 1:length(tvec)
    xtrue(:,i) = ppval(x,tvec(i));
    epsdiff(i)=abs(epshat(i)-Trueeps(i));
    xdiff(i) = norm(xtrue(:,i)-xhat(:,i));
end
%xdiff = zeros(size(tvec));


figure(9)
plot(tvec,xdiff,'b--', tvec, epsdiff, 'r-.','LineWidth',2)
yline(0);
xline(Ldouble);
xlim([0,Tstar]);
legend('$|\hat{x}(t)-x(t)|$','$|\hat{\epsilon}(t)-\epsilon|$','Interpreter','latex','FontSize',11);
ax = gca;
ax.FontSize=20;

figure(10)
plot(tvec,xdiff,'b--', tvec, epsdiff, 'r-.','LineWidth',2)
yline(0);
xlim([Ldouble,Tstar]);
legend('$|\hat{x}(t)-x(t)|$','$|\hat{\epsilon}(t)-\epsilon|$','Interpreter','latex','FontSize',11);
ax = gca;
ax.FontSize=20;

% figure(3)
% plot(t,x1plot,'b-','LineWidth',1.5);
% hold on 
% plot(tvec, xhat1, 'b-.','LineWidth', 2);
% hold on 
% plot(t, x2plot,'r-','LineWidth',1.5);
% hold on 
% plot(tvec,xhat2,'r--','LineWidth',2);
% hold on 
% xline(Ldouble)
% legend('$x_1(t)$', '$\hat{x}_1(t)$', '$x_2(t)$', '$\hat{x}_2(t)$', '2L','Interpreter', 'latex', 'FontSize', 11)
% xlim([0,Tstar])
% ax = gca;
% ax.FontSize=20;
% axes('Position',[.2,.7,.2,0.7])
% box on 
% plot(t,x1plot,'b-',tvec,xhat1,'b-','LineWidth',2);
% hold on 
% plot(t,x2plot,'r-',tvec,xhat2,'r--','LineWidth',2);
% xlim([Ldouble,Tstar])

end