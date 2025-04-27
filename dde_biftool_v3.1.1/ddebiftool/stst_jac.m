function [J,res]=stst_jac(funcs,x,par,free_par)
%% residual and jacobian for equilibrium problem
% function [J,res]=stst_jac(funcs,x,par,free_par)
% INPUT:
%   funcs problem functions
%	x current solution guess in R^n
%	par current parameter values
%	free_par free parameter numbers
% OUTPUT: 
%	J jacobian in R^(n+s x n+p)
%	res residual in R^(n+s x 1)

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
sys_tau=funcs.sys_tau;
sys_rhs=funcs.sys_rhs;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

n=length(x); % system dimension 

if funcs.tp_del==0,
    m=length(sys_tau()); % number of delays
else
    m=sys_ntau();
end

xx=x(:,ones(m+1,1));

res=sys_rhs(xx,par);

J=zeros(n,n+length(free_par));

for i=0:m
    J(1:n,1:n)=J(1:n,1:n)+sys_deri(xx,par,i,[],[]);
end;
for j=1:length(free_par)
    J(1:n,n+j)=J(1:n,n+j)+sys_deri(xx,par,[],free_par(j),[]);
end
end
