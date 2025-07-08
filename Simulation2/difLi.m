function dLi = difLi(t,y,yL,tau)
rho1 = [0;1];
A = [0,1;0,0];
C = [1,0];
y1 = y(1:2);
y2 = y(3:end);


dy1 = A*y1+rho1;
dy2 = -transpose(A)*y2+transpose(C)*C*expm(tau*A)*yL(1:2);
dLi = [dy1;dy2];
end