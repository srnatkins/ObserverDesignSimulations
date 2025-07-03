function [Beta11,Beta12,Beta13] =interpolatebeta1(tint,k,tau,h)
[t,beta11,beta12,beta13]= beta1solve(tint,k,tau,h);
Beta11 = spline(t,beta11);
Beta12 = spline(t,beta12);
Beta13 = spline(t,beta13);
end