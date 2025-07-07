function x = interpolate_x(t,eps,x0)
[tnew,x]= xsolve(t,eps,x0);
x=transpose(x);
x = spline(tnew,x);
end