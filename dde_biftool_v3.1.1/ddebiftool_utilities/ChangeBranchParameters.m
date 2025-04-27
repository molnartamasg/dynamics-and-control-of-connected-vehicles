function [newbranch,suc]=ChangeBranchParameters(funcs,branch,ind,varargin)
%% Change parameters for continuation along branch and create first two points
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
default={'contpar',[],'dir',[],'step',0.01,'correc',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
newbranch=replace_branch_pars(branch,options.contpar,pass_on);
pini=branch.point(ind);
if isempty(options.dir)
    options.dir=newbranch.parameter.free(1);
end
[newbranch,suc]=correct_ini(funcs,newbranch,pini,...
    options.dir,options.step,options.correc);
end
