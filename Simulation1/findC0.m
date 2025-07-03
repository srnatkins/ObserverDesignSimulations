function C0 = findC0(C,tau)
A =[0 1 0 0 0; -1 -1 1 0 0; 0 0 0 1 0; 1 0 -1 -1 1; 0 0 0 -1 -1];
M= expm(A*tau);
 if det(M)==0
    disp('Not invertible')
 else
 C0 = C*M;
 end
end