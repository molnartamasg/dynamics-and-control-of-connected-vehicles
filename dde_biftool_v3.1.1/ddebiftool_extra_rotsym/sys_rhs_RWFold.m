function y=sys_rhs_RWFold(x,p,rho,orig_rhs,dim,hbif)
%% rhs of extended DDE for fold of relative equilibria
%
% x extended state (orig state and null vector)
% p user parameters
% rho rotation frequency
% orig_rhs user r.h.s
% dim original system dimension
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
x0=x(1:dim,:);
v=x(dim+1:end,:);
y0=orig_rhs(x0,p);
dfx=@(x0,dev,ind)app_dir_deriv(@(x)orig_rhs(x,p),x0,dev,ind,hbif);
%% add partial derivatives of all terms
y1=0*y0;
for i=1:size(x0,2)
    y1=y1+dfx(x0,v(:,i),i);
end
dfom=app_dir_deriv(@(p)orig_rhs(x0,p),p,1,length(p),hbif);
y1=y1+dfom*rho;
y=[y0;y1];
end
