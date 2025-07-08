function [Tstar,tint,L]=getTstar(tau,tau1,h,gap)
% We want to ensure time interval time interval tint=[0,T_*] is set to
% where T_*>2L where L=max_i(tau_i)+h+tau as specified in Theorem 2
% moreover we want time interval Tstar-2L=gap so that time interval
% [2L,T_*) is large enough to study

L= tau1+h+tau;
Tstar=2*L+gap;
tint=[0,Tstar];

end
