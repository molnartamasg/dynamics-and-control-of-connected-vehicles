function [r,J]=sys_cond_MWFold(point,usercond,...
    dim,beta_ind,period_ind,rho_ind,relations,get_comp)
%% constraints used for extended DDE in fold continuation for relative periodic orbits
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%% user conditions
Jtemplate=p_axpy(0,point,[]);
userpoint=get_comp(point,'solution');
nuserpar=length(userpoint.parameter);
[userres,userJ]=usercond(userpoint);
userJext=repmat(Jtemplate,length(userJ),1);
%% append artificial parameters and components
for i=1:length(userJ)
    userJext(i).parameter(1:nuserpar)=userJ(i).parameter;
    userJext(i).profile(1:dim,:)=userJ(i).profile;
    userJext(i).period=userJ(i).period;
end
%% duplicate rotational phase condition (last user condition)
%residual is set to 0 but Jacobian pretends it is applied to the nullvector
userres(end+1)=0;
userJext(end+1)=userJext(end);
userJext(end).profile=userJext(end).profile([dim+1:end,1:dim],:);
%% obtain condition <v,v>+beta^2+rho^2-1=0
vpoint=p_axpy(0,point,[]);
vpoint.profile(dim+1:end,:)=point.profile(dim+1:end,:);
vpoint.parameter=point.parameter;
[vtv,vtvJ]=p_dot(vpoint,vpoint,'free_par_ind',[beta_ind,rho_ind]);
vtvJ=p_axpy(2,vtvJ,[]);
vtvres=vtv-1;
%% obtain condition <x',v>=0
xdot=point;
xdot.profile=[zeros(dim,size(xdot.profile,2));xdot.profile(1:dim,:)];
[xdtvres,xdtvJv]=p_dot(point,xdot,'derivatives',[0,1]);
[xdtvresdum,xdtvJx]=p_dot(xdot,point,'derivatives',[1,0]); %#ok<ASGLU>
xdtvJ=xdtvJv;
xdtvJ.profile(1:dim,:)=xdtvJx.profile(dim+1:2*dim,:);
%% constrain parameter(period_ind)-period=0;
pres=point.parameter(period_ind)-point.period;
periodJ=p_axpy(0,point,[]);
periodJ.parameter(period_ind)=1;
periodJ.period=-1;
%% fix relations between extended delays and original delays
taures=relations*point.parameter(:);
tauJ=repmat(Jtemplate,size(relations,1),1);
for j=1:size(relations,1)
    tauJ(j).parameter=relations(j,:);
end
%% assemble residuals and Jacobians
r=[userres(:);vtvres;xdtvres;pres;taures];
J=[userJext(:);vtvJ;xdtvJ;periodJ;tauJ];
end
