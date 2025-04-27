function f = cusp_rhs(xx,par)
%% right-hand side of cusp demo (see cusp_demo.m for equations)
%
% (c) DDE-BIFTOOL v. 3.1.1(121), 02/09/2015
%
%par = [q11,q12,q21,e1,e2]

% f(1,1) = -xx(1,1)+par(1)/(1+exp(-4*xx(1,2)))-par(2)*xx(2,2)+par(4);
% f(2,1) = -xx(2,1)+par(3)/(1+exp(-4*xx(1,2)))+par(5);

alpha = @(u) 1/(1+exp(-4*u))-1/2;

f(1,1) = -xx(1,1)+par(1)*alpha(xx(1,2))-par(2)*xx(2,2)+par(4);
f(2,1) = -xx(2,1)+par(3)*alpha(xx(1,2))+par(5);
end
