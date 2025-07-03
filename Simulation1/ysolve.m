function y = ysolve(tint,eps,k,x0,tau)
C = [1, 0, 0, 0, 0];
C0 = findC0(C,tau);
x = interpolate_x(tint,eps,k,x0);
t = 0:1.e-3:tint(end);
y = zeros(size(t));
db = @(t) 3*sin(t);
for i = 1:length(t)
    if t(i)<tau
        y(i)=C0*x0+db(t(i));
    else
        y(i) = C0*ppval(x,t(i)-tau)+db(t(i));
    end
end
y = spline(t,y);

end