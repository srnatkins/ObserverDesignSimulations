function main
A=[0 1; 0 0];
C = [1 0];
h=3.8;
tau1=0;
tau=0.125;
gap=11.9+20.25;
eps =1;
u0=[1;0;1;0];
[E,E1]=findE(h);
[Tstar,tint,L]=getTstar(tau,tau1,h,gap);
Ldouble = 2*L;
fprintf('L=h+tau1+tau=%f.\n',L);
fprintf('T_star is %f.\n',Tstar);
fprintf('2L is %f.\n', Ldouble);
fprintf('T_star>2L since T*-2L=%f.\n',gap);
[betastar, Lstar, betabarnew,smallkd,largekd] = NewAssumption5_1(tint,tau,h,tau1);
kd= smallkd;
[tvec,xhat,epshat] = getxhatandepshat(tint,eps,u0,kd,Lstar, tau,tau1,h);
Trueeps=eps*ones(size(tvec));
[xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd);
xtrue = zeros(2,length(tvec));
epsdiff = zeros(size(tvec));
xdiff = zeros(size(tvec));
for i = 1:length(tvec)
    xtrue(:,i) = ppval(x,tvec(i));
    epsdiff(i)=abs(epshat(i)-Trueeps(i));
    xdiff(i) = norm(xtrue(:,i)-xhat(:,i));
end
xdiffsmall = xdiff;
epsdiffsmall = epsdiff;
tvecsmall = tvec;


kd=largekd;
[tvec,xhat,epshat] = getxhatandepshat(tint,eps,u0,kd,Lstar, tau,tau1,h);
Trueeps=eps*ones(size(tvec));
[xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd);
xtrue = zeros(2,length(tvec));
epsdiff = zeros(size(tvec));
xdiff = zeros(size(tvec));
for i = 1:length(tvec)
    xtrue(:,i) = ppval(x,tvec(i));
    epsdiff(i)=abs(epshat(i)-Trueeps(i));
    xdiff(i) = norm(xtrue(:,i)-xhat(:,i));
end
xdiffbig = xdiff;
epsdiffbig = epsdiff;
tvecbig = tvec;

figure(1)
plot(tvecsmall,xdiffsmall,'b-', tvecbig, xdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([0,Tstar]);
legend('$|\hat{x}(t)-x(t)|$ with $k_{\Delta}=0.121364$', '$|\hat{x}(t)-x(t)|$ with $k_{\Delta}=0.236977$','$2L$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

figure(2)
plot(tvecsmall,xdiffsmall,'b-', tvecbig, xdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([Ldouble,Tstar]);
legend('$k_{\Delta}=0.121364$', '$k_{\Delta}=0.236977$','Interpreter','latex','FontSize',11)
ylabel('$|\hat{x}(t)-x(t)|','Interpreter','latex','FontSize',20)
ax=gca;
ax.FontSize=20;

figure(3)
plot(tvecsmall,epsdiffsmall,'b-', tvecbig, epsdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([0,Tstar]);
legend('$|\hat{\epsilon}(t)-\epsilon|$ with $k_{\Delta}=0.121364$', '$|\hat{\epsilon}(t)-\epsilon|$ with $k_{\Delta}=0.236977$','$2L$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

figure(4)
plot(tvecsmall,epsdiffsmall,'b-', tvecbig, epsdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([Ldouble,Tstar]);
legend('$|\hat{\epsilon}(t)-\epsilon|$ with $k_{\Delta}=0.121364$', '$|\hat{\epsilon}(t)-\epsilon|$ with $k_{\Delta}=0.236977$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;


figure(5)
plot(tvecsmall,xdiffsmall,'b-', tvecbig, xdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([0,Tstar]);
ylabel('$|\hat{x}(t)-x(t)|$','Interpreter','latex','FontSize',20)
legend('$k_{\Delta}=0.121364$', '$k_{\Delta}=0.236977$','$2L$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

figure(6)
plot(tvecsmall,xdiffsmall,'b-', tvecbig, xdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([Ldouble,Tstar]);
ylabel('$|\hat{x}(t)-x(t)|$','Interpreter','latex','FontSize',20)
legend('$k_{\Delta}=0.121364$', '$k_{\Delta}=0.236977$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

figure(7)
plot(tvecsmall,epsdiffsmall,'b-', tvecbig, epsdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([0,Tstar]);
ylabel('$|\hat{\epsilon}(t)-\epsilon|$', 'Interpreter','latex','FontSize',20)
legend('$k_{\Delta}=0.121364$', '$k_{\Delta}=0.236977$','$2L$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

figure(8)
plot(tvecsmall,epsdiffsmall,'b-', tvecbig, epsdiffbig,'r--','LineWidth',2)
xline(Ldouble);
xlim([Ldouble,Tstar]);
ylabel('$|\hat{\epsilon}(t)-\epsilon|$', 'Interpreter','latex','FontSize',20)
legend('$k_{\Delta}=0.121364$', '$k_{\Delta}=0.236977$','Interpreter','latex','FontSize',11)
ax=gca;
ax.FontSize=20;

end
