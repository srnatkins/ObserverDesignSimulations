function [L11, L21] = interpolateL(t,tau)
[t1, L11, L21,sol1] = Lsolve(t,tau);
L11 = spline(t1,L11);
L21 = spline(t1,L21);
end