function x = interpolate_x(t,eps,k,x0)
[tnew,x]= xsolve(t,eps,k,x0);
x=transpose(x);
x = spline(tnew,x);
end