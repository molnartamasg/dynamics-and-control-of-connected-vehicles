%% Branch off at period doubling (branch psolbranch, point number ind)
%%
function [per2,suc]=DoublePsol(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |psolbranch|: branch of |'psol'| periodic orbits from which one wants to branch off
% * |ind|: index in |'point'| field of |psolbranch| which is close to
% period doubling where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'radius'|: initial deviation along period-doubling eigenvector
%
% All other name-value pairs are passed on to output branch.
%% Outputs
%
% * |per|: branch of periodic orbits with desired settings, double period,
% and two initial corrected points
% * |suc|: flag indicating success
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%%
default={'radius',0.01};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch per2 of periodic solutions starting from an
% approximate period doubling num on a branch br of periodic orbits
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
    return
end
[psol2,suc]=p_correc(funcs,deg_psol2,per2.parameter.free,step_cond2,...
    per2.method.point,1,deg_psol2);
if suc==0
    per2=[];
    return
end
per2.point=psol1;
per2.point(2)=psol2;
end
