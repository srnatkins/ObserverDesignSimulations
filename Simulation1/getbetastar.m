function [betastar,t,beta11,beta12,beta13] = getbetastar(tint,k,tau,h,tau1,tau2,tau3)
[t,Beta31,Beta32,Beta33,beta11,beta12,beta13] = beta3solve(tint,k,tau,h);
dt = 0.05;
istar1 = round(tau1/dt)+1
istar2 = round(tau2/dt)+1
istar3 = round(tau3/dt)+1
% for i = 1:length(t)
%     if t(i) == tau1
%         istar1 = i;
%         disp('istar1 found')
%     elseif t(i) == tau2
%         istar2 = i;
%         disp('istar2 found')
%     elseif t(i) == tau3
%         istar3 =i;
%         disp('istar3 found')
% %     else
% %         istar3 = 21;
% %         disp(t(istar3))
%     end 
% end
beta31tau1 = Beta31(istar1);
beta31tau2 = Beta31(istar2);
beta31tau3 = Beta31(istar3);
beta32tau1 = Beta32(istar1);
beta32tau2 = Beta32(istar2);
beta32tau3 = Beta32(istar3);
beta33tau1 = Beta33(istar1);
beta33tau2 = Beta33(istar2);
beta33tau3 = Beta33(istar3);
beta3tau1 = [beta31tau1,beta32tau1,beta33tau1];
beta3tau2 = [beta31tau2,beta32tau2,beta33tau2];
beta3tau3 = [beta31tau3,beta32tau3,beta33tau3];
betastar = [beta3tau1;beta3tau2;beta3tau3];
end
