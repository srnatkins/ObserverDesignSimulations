function C0 = findC0(C,tau)
A =[0 1; 0 0];
M= expm(A*tau);
 if det(M)==0
    disp('Not invertible')
 else
 %Minv = inv(M);
 %C0 = C*Minv;
 %c0 = transpose(M)\transpose(C);
 %C0 = transpose(c0);
 C0 = C*M;
 end
end