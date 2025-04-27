function [r,J]=sys_cond_TorusBif(point,usercond,dim,period_ind,get_comp)
%% additional conditions for torus and period doubling bifurcation
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
Jtemplate=p_axpy(0,point,[]);
%% user conditions
userpoint=get_comp(point,'solution');
[userres,userJ]=usercond(userpoint);
nuserpar=length(userpoint.parameter);
userJext=repmat(Jtemplate,length(userJ),1);
%% append artificial parameters and components to user conditions
for i=1:length(userJ)
    userJext(i).parameter(1:nuserpar)=userJ(i).parameter;
    userJext(i).profile(1:dim,:)=userJ(i).profile;
    userJext(i).period=userJ.period;
end
%% keep Floquet mode at unit length
uvpoint=Jtemplate;
uvpoint.profile(dim+1:end,:)=point.profile(dim+1:end,:);
[utuvtv,utuvtvJ]=p_dot(uvpoint,uvpoint);
utuvtvJ=p_axpy(2,utuvtvJ,[]);
utuvtvres=utuvtv-1;
%% ensure that real and imaginary part of Floquet mode are orthogonal
vupoint=uvpoint;
vupoint.profile(dim+1:3*dim,:)=uvpoint.profile([2*dim+1:3*dim,dim+1:2*dim],:);
[utv,utvJ]=p_dot(uvpoint,vupoint);
utvres=utv/2;
%% constrain parameter(period_ind)-period=0;
pres=point.parameter(period_ind)-point.period;
periodJ=Jtemplate;
periodJ.parameter(period_ind)=1;
periodJ.period=-1;%% keep 
r=[userres;utuvtvres;utvres;pres];
J=[userJext;utuvtvJ;utvJ;periodJ];
end
