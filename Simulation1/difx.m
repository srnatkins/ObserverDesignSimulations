function dx = difx(t,x,eps,k)
    %define basis functions
    lam = [1,2,3];
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
    rho2 = @(t)  [f( (t-lam(2))^2/( 2*w(2)^2 ) );f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) ); f( (t-lam(2))^2/( 2*w(2)^2 ) )];
    rho3 = @(t) [f( (t-lam(3))^2/( 2*w(3)^2 ) );f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) ); f( (t-lam(3))^2/( 2*w(3)^2 ) )];
    da = @(t) 3*[1*sin(t); 2*sin(2*t); 3*sin(3*t); 4*sin(4*t);5*sin(5*t)];
    A = [0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
    dx = A*x+eps(1)*rho1(t)+eps(2)*rho2(t)+eps(3)*rho3(t)+da(t);
end
