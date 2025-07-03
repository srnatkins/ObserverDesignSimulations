function [gam1,gam2,gam3] = gamsolve(k,tau)
%defining A and C
A = [0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
C = [1 0 0 0 0];
%get rho1,rho2, rho3
lam =[1,2,3];
w =[1,2,3];
if k == 1 %Guassian form
    f = @(t) exp(-t);
elseif k == 2 %multiquadratic
    f = @(t) sqrt(1+t);
elseif k == 3 %inverse multiquadratic
    f = @(t) (1+t)^(-1/2);
elseif k == 4 %cauchy form
    f = @(t) (1+t)^(-1);
end
rho1 = @(t) [f( (t-lam(1))^2/( 2*w(1)^2 ) );f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) ); f( (t-lam(1))^2/( 2*w(1)^2 ) )];
rho2 = @(t) [f( (t-lam(2))^2/( 2*w(2)^2 ) );f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) )];
rho3 = @(t) [f( (t-lam(3))^2/( 2*w(3)^2 ) );f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) )];
%find gam1,gam2,gam3
fun1 = @(s) expm(A*(-s))*rho1(s);
fun2 = @(s) expm(A*(-s))*rho2(s);
fun3 = @(s) expm(A*(-s))*rho3(s);
gam1 = @(t) -C*expm(A*t)*integral(fun1,t-tau,t,'ArrayValued',true);
gam2 = @(t) -C*expm(A*t)*integral(fun2,t-tau,t,'ArrayValued',true);
gam3 = @(t) -C*expm(A*t)*integral(fun3,t-tau,t,'ArrayValued',true);

end
