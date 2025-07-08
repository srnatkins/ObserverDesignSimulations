function [E,E1] = findE(h)
A =[0 1; 0 0];
C = [1 0];

fun =@(s) expm(transpose(A)*s)*transpose(C)*C*expm(A*s);
E = integral(fun,-h,0,'ArrayValued',true);

E1 = inv(E);
end