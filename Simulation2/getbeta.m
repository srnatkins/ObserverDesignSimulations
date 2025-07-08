function betastar = getbeta(t,tau,tau1,h)
[L11, L21] = interpolateL(t,tau);
A =[0 1; 0 0];
C = [1 0];
[E,E1] = findE(h);  %E1 is the inverse of E
if tau1> tau+h  | tau1==h+tau %if tau1 satisfies condition described in remark 3.7 of manuscript
    betastar = C*(expm(tau*A)*ppval(L11,tau1-tau))-C*E1*(ppval(L21,tau1)-expm(-h*transpose(A))*ppval(L21,tau1-h)); %beta*=beta**=beta31(tau1)
else
    time = 0:0.1:t(end);
    betastarvec = zeros(size(time));
    for i = 1: length(time)
        betastarvec(i) = C*(expm(tau*A)*ppval(L11,time(i)-tau1-tau))-C*E1*(ppval(L21,time(i)-tau1)-expm(-h*transpose(A))*ppval(L21,time(i)-tau1-h));%betastarvec(i)=beta31(t(i)-tau1).  
    end
    % When t(i)>h+tau1+tau, betastarvec(i) is constant.
    % we set betastar = betastarvec(end)=beta31(T_*-tau1)
   betastar=betastarvec(end);  
end





end