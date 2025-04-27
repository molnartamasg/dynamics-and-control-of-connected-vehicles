function y=sys_rhs_POfold(x,p,beta,period,tau,sys_rhs,dim,xtau_ind,sys_deri)
%% rhs of extended DDE for fold of periodic orbits
% tau contains delays ([0,p(sys_tau())]), 
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%

xall=x(1:dim,:,:);
x0=xall(:,xtau_ind(1,:),:);
v=x(dim+1:end,xtau_ind(1,:),:);
y0=sys_rhs(x0,p);
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)sys_rhs(x,p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(sys_deri(x0,p,ind-1,[],[]),dev,'*');
end
%% partial derivative wrt non-delayed term
y1=beta/period*y0+df(x0,v(:,1,:),1);
%% add partial derivatives of all delayed terms
for i=2:size(x0,2)
    deviation=v(:,i,:)+sys_rhs(xall(:,xtau_ind(i,:),:),p)*tau(i)*beta/period;
    y1=y1+df(x0,deviation,i);
end
y=cat(1,y0,y1);
end
