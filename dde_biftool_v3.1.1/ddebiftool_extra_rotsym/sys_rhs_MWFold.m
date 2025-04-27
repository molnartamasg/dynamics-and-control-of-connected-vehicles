function y=sys_rhs_MWFold(x,p,beta,period,rho,tau,sys_rhs,dim,xtau_ind,hbif)
%% rhs of extended DDE for fold of modulated waves
%
% tau contains delays ([0,p(sys_tau())]), 
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% beta is derivative wrt period
% rho is rotation frequency
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%

xall=x(1:dim,:,:);
x0=xall(:,xtau_ind(1,:),:);
v=x(dim+1:end,xtau_ind(1,:),:);
y0=sys_rhs(x0,p);
df=@(x0,dev,ind)app_dir_deriv(@(x)sys_rhs(x,p),x0,dev,ind,hbif);
%% partial derivative wrt non-delayed term
y1=beta/period*y0+df(x0,v(:,1,:),1);
%% add partial derivatives of all delayed terms
for i=2:size(x0,2)
    deviation=v(:,i,:)+...
        sys_rhs(xall(:,xtau_ind(i,:),:),p)*tau(i)*beta/period;
    y1=y1+df(x0,deviation,i);
end
dfom=app_dir_deriv(@(p)sys_rhs(x0,p),p,1,length(p),hbif);
y1=y1+dfom*rho;
y=cat(1,y0,y1);
end
