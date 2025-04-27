function y=sys_rhs_TorusBif(x,p,omega,period,tau,sys_rhs,dim,sys_deri)
%% r.h.s. of extended DDE for torus and period doubling bifurcation
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
x0=x(1:dim,:,:);
u=x(dim+1:2*dim,:,:);
v=x(2*dim+1:end,:,:);
y0=sys_rhs(x0,p);
yu=pi*omega/period*v(:,1,:);
yv=-pi*omega/period*u(:,1,:);
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)sys_rhs(x,p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(sys_deri(x0,p,ind-1,[],[]),dev,'*');
end
for i=1:size(x0,2)
    %% add partial derivatives of all delayed terms (incl delay zero)
    c=cos(pi*omega/period*tau(i));
    s=sin(pi*omega./period*tau(i));
    yu=yu+df(x0, c*u(:,i,:)+s*v(:,i,:),i);
    yv=yv+df(x0,-s*u(:,i,:)+c*v(:,i,:),i);
end
y=cat(1,y0,yu,yv);
end
