function y = ysolve(tint,eps,x0,tau)
C=[1,0];
C0 = findC0(C,tau);
x = interpolate_x(tint,eps,x0);
t = 0:1.e-3:tint(end);
y = zeros(size(t));
for i = 1:length(t)
    if t(i)<tau
        y(i)=C0*x0;
    else
        y(i) = C0*ppval(x,t(i)-tau);
    end
end
y = spline(t,y);
end