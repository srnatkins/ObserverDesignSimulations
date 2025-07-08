function [Dsharp,Ddoublestar]= getDs(tint,eps,u0,kd,tau,tau1,h)
 A=[0 1; 0 0];
 C=[1,0];
 Csharp = transpose(C)*C;
 [E,E1] = findE(h);  %E1 is the inverse of E
 [Deld, Dd] = getDeldDd(tint,eps,u0,kd,tau,h);
 func = @(m) expm(-A*m)*Deld(m);
 L4= @(s,t) expm(A*s)*integral(func,s-tau,t,'ArrayValued',true);
 func2 = @(s,t) E1*expm(transpose(A)*(s-t))*Csharp*L4(s,t);
 %time = tint(1): 0.01: tint(end);
 Dsharp = @(t) integral(@(s)func2(s,t),t-h,t,'ArrayValued',true);
 Dstar = @(t) C*Dsharp(t)+Dd(t);
 Ddoublestar= @(t) Dstar(t-tau1);
end