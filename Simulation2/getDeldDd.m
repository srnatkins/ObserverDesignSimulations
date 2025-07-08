function [Deld, Dd] = getDeldDd(tint,eps,u0,kd,tau,h)
%based upon diffu \Delta(\xi(t),t)=[kdsin(\xi_1(t)),0]
% we had manually input \Delta and \Delta_d into diffu in order to
% numerically compute
% solving system [\hat{\xi},x]. Now we want to take the numerical solution
% obtained and extract what \Delta_d is and \mathcal{D}_d which are described in
% equation (32)

%output 
%       Deld -> \Delta_d(x(t),t)
%       Dd   -> \mathcal{D}_d(x_t,t)
[E,E1] = findE(h);
C=[1,0];
A = [0 1; 0 0];
[xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd);
Del1 = @(t) [kd*sin(ppval(x1,t)+ppval(xihat1,t));0];
Del2 = @(t) [kd*sin(ppval(xihat1,t));0];
Deld = @(t) Del1(t)-Del2(t);
%C0 = findC0(C,tau);
func1 = @(ell) expm(-A*ell)*Deld(ell);
Dd = @(t) C*expm(A*t)*integral(func1,t-tau,t,'ArrayValued',true);
% time = tint(1):0.1:tint(end);
% dd = zeros(size(time));
% for i = 1:length(time)
%     func = @(l) -C0*expm(A*(time(i)-l-tau))*Deld(l);
%     dd(i) = integral(func,time(i)-tau,time(i),'ArrayValued',true);
% end

%Dd = spline(time,dd);


% [t,u,xihat,x] = usolve(tint,eps,u0,kd);
% xihat = transpose(xihat);
% x = transpose(x);
% deld=zeros(size(x));
% deld(1,:) = kd*sin(x(1,:)+xihat(1,:))-kd*sin(xihat(1,:));
% Deld = spline(t,deld);
% C0 = findC0(C,tau);
% [xihat, x] = interpolate_xihatandx(tint,eps,u0,kd);
% tnew = tint(1): 0.1: tint(end);
% n= length(tnew);
% dd = zeros(2,n)
% for i = 1:n
%     fun1 = @(l) expm(A*(tnew(i)-l-tau))*ppval(Deld())
% end
end
