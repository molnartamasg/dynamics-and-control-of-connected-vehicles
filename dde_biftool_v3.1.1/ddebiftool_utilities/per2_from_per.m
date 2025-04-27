function [per2,suc]=per2_from_per(funcs,psolbranch,ind,varargin)
%% branch off at period doubling (branch psolbranch, point number ind)
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%%
default={'radius',0.01};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per of periodic solutions starting from an
% approximate Hopf point num on a branch br of steady states
point=psolbranch.point(ind);
intervals=(size(point.profile,2)-1)/point.degree;
[deg_psol1,step_cond1]=p_topsol(funcs,point,options.radius,...
    [],intervals);
[deg_psol2,step_cond2]=p_topsol(funcs,point,options.radius*2,...
    [],intervals);
per2=replace_branch_pars(psolbranch,[],pass_on);
[psol1,suc]=p_correc(funcs,deg_psol1,per2.parameter.free,step_cond1,...
    per2.method.point,1,deg_psol1);
if suc==0
    per2=[];
    return;
end;
[psol2,suc]=p_correc(funcs,deg_psol2,per2.parameter.free,step_cond2,...
    per2.method.point,1,deg_psol2);
if suc==0
    per2=[];
    return;
end;
per2.point=psol1;
per2.point(2)=psol2;
end
