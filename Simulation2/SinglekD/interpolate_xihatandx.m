function [xihat, x,xihat1,x1] = interpolate_xihatandx(tint,eps,u0,kd)
%returns interpolated x and \hat\xi via spline interpolation
[tnew,u,xihat,x] = usolve(tint,eps,u0,kd);
xinodes = transpose(xihat);
xnodes = transpose(x);
%interpolate xihat and x to where xihat(t)=xihat(0) and x(t)=x(0) for t<0
tneg= -tint(end):0.1:-.1;
xihist1 = u0(1)*ones(size(tneg));
xihist2 = u0(2)*ones(size(tneg));
xi_hist =[xihist1;xihist2];
xhist1 = u0(3)*ones(size(tneg));
xhist2 = u0(4)*ones(size(tneg));
x_hist = [xhist1;xhist2];
tnew = tnew';
tfull = [tneg,tnew];
xi_nodes = [xi_hist,xinodes];
x_nodes = [x_hist,xnodes];
xihat1 = spline(tfull,xi_nodes(1,:));
x1 = spline(tfull,x_nodes(1,:));
xihat = spline(tfull,xi_nodes);
x = spline(tfull,x_nodes);



end
