function [r,J]=sys_cond_RWFold(p,orig_cond,dim,ind_rho)
%% constraints used for extended DDE in fold continuation of relative equilibria
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
p_orig=p;
p_orig.x=p.x(1:dim);
p_orig.parameter=p_orig.parameter(1:ind_rho-1);
[r_orig,J_orig]=orig_cond(p_orig);
nuser=length(r_orig);
for i=1:nuser
    J_orig(i).x=[J_orig(i).x;zeros(dim,1)];
    J_orig(i).parameter=[J_orig(i).parameter,0];
end
rphas=0;
Jphas=J_orig(end);
Jphas.x=[zeros(dim,1);Jphas.x(1:dim)];
rnorm=sum(p.x(dim+1:end).^2)+p.parameter(ind_rho)^2-1;
Jnorm=p_axpy(0,p,[]);
Jnorm.x(dim+1:end)=2*p.x(dim+1:end);
Jnorm.parameter(ind_rho)=2*p.parameter(ind_rho);
r=[r_orig(:);rphas;rnorm];
J=[J_orig(:);Jphas;Jnorm];
end
