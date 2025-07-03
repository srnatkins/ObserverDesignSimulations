function [t,Beta31,Beta32,Beta33,beta11,beta12,beta13] = beta3solve(tint,k,tau,h)
C = [1,0,0,0,0];
%find beta1i
[t,beta11,beta12,beta13]= beta1solve(tint,k,tau,h);
%[Beta11,Beta12,Beta13] =interpolatebeta1(tint,k,tau,h);
[gam1,gam2,gam3] = gamsolve(k,tau);
Beta31 = zeros(size(t));
Beta32 = zeros(size(t));
Beta33 = zeros(size(t));
for i =1:length(t)
    Beta31(i) =  C*beta11(:,i)+gam1(t(i));
    Beta32(i) =  C*beta12(:,i)+gam2(t(i));
    Beta33(i) =  C*beta13(:,i)+gam3(t(i));
end
end
