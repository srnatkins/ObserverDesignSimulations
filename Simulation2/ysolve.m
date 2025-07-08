function y = ysolve(tint,eps,u0,kd,tau)
C = [1 0];
C0 = findC0(C,tau);
[xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd);
t = tint(1): 1.e-3:tint(end);
y = zeros(size(t));
for i = 1:length(t)
    if t(i)<tau
        y(i)=C0*u0(3:4);
    else
        y(i)=C0*ppval(x,t(i)-tau);
    end
end
y = spline(t,y);

end