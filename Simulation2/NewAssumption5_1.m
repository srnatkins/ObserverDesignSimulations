function [betastar, Lstar, betabarnew,smallkd,largekd] = NewAssumption2_3(tint,tau,h,tau1)
%Function returns betastar,Lstar, betabarnew, and smallkd and largekd for simulation 2
%This is needed in order to generate the comparison plots between two
%Lipschitz functions \Delta(\xi,t) with Lipschitz constants being smallkd
%and largekd respectively (where those constants specifically satisfy
%Assumption 5.1. 
%INPUTS: 
%       tint        time interval [0,T_*]
%       tau         delay
%       h    
%       tau1
%%NOTE it is imporitant that inputs are chosen to where
% T_^*>h+tau1+tau hold as this is one condition needed for Assumption 2.2 
%OUTPUTS:
%       betastar    constructed based upon assumption 2
%                   where beta3 is constructed via Remark 12  
%       Lstar       the left inverse of betastar           
%       betabarnew  constructed via equation (67)
%       smallkd     smallkd = 1/(2betabarnew) 
%       largekd     largekd = 1/(betabarnew+0.1)     

%Based upon sim. 2 computations, betastar is a constant function so
%betastar is returned as number. Additionally, Lstar=1/betastar is a number 
%which exists so long as betastar is nonzero.  This is along with 
% T_*>h+tau1 will make it to where Assumption 2.2 to be satisfied

%kd>0 chosen such that kd*betabarnew<1 which is a condition for 
% Assumptioon 5.1. 

A=[0 1; 0 0];
C=[1,0];
Csharp = transpose(C)*C;
[E,E1] = findE(h);  %E1 is the inverse of E
betastar = getbeta(tint,tau,tau1,h);
betacheck = 1/12*h^2;
fprintf('beta_star is %f.\n',betastar)
fprintf('h^2/12 is %f.\n',betacheck)
Lstar= 1/betastar;
[t,beta1] =beta1solve(tint,tau,h);
fprintf('L* is %f.\n', Lstar);

M1 = beta1.*Lstar;
fun1= @(l)norm(M1*C*expm(-A*l));
fun2 = @(s,m) norm((eye(2)-M1*C)*E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
term1 = integral(fun1,-tau,0,'ArrayValued',true);

term2 = integral(@(s) integral(@(m)fun2(s,m),s-tau,0,'ArrayValued',true),-h,0,'ArrayValued',true);
betabarnew = term1+term2;
Oneoverbetabar=1/betabarnew;
fprintf('betabar_new is %f.\n', betabarnew);
fprintf('1/betabar_new is %f.\n', Oneoverbetabar);

smallkd=1/(2*betabarnew);
assump3=smallkd*betabarnew;
fprintf('Setting k_Delta as 1/(2*betabar_new) yields %f. \n',smallkd);
fprintf('Thus k_Delta*betabar_new is %f. \n', assump3);

largekd =1/(betabarnew+0.1);
assump3 = largekd*betabarnew;
fprintf('Setting k_Delta as 1/(betabar_new+0.1) yields %f. \n',largekd);
fprintf('Thus k_Delta*betabar_new is %f. \n', assump3);
% func1 = @(l)norm(M1*C*expm(A*l));
% func2 = @(s,m)norm(M1*C*E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
% func3 = @(s,m)norm(E1*expm(transpose(A)*s)*Csharp*expm(A*(s-m)));
% Term1 = integral(func1,0,tau,'ArrayValued',true)
% Term2 = integral(@(s) integral(@(m)func2(s,m),s-tau,0,'ArrayValued',true),-h,0,'ArrayValued',true)
% Term3 = integral(@(s) integral(@(m)func3(s,m),s-tau,0,'ArrayValued',true),-h,0,'ArrayValued',true)

% betabar=Term1+Term2+Term3
% kd=1/(betabar+0.1)
%kdcheck=1/betabar -0.1

end