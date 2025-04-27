function branch=replace_branch_pars(branch,contpar,pass_on)
%% pass on optional arguments to branch structure
%
% branch=replace_branch_pars(branch,contpar,pass_on)
% 
% input
% branch: given branch structure to be amended
% contpar: continuation parameter, prepended to free parameters already
%          present in branch (if length==1) or replacing
%          branch.parameter.free
% pass_on: cell array containing name-value pairs for fields of
%          branch.method.continuation, branch.method.point, branch.method.stability and
%          branch.parameter to be replaced
%
% output: branch with amended fields
%
% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
%%
branch.method.continuation=dde_set_options(branch.method.continuation,pass_on,'pass_on');
branch.method.point=dde_set_options(branch.method.point,pass_on,'pass_on');
if isfield(branch.method,'stability')
    branch.method.stability=dde_set_options(branch.method.stability,pass_on,'pass_on');
end
if isempty(contpar)
    % if no continuation parameters are given use free parameters of branch
    branch.parameter.free=branch.parameter.free;
else
    branch.parameter.free=contpar(:)';
end
branch.parameter=dde_set_options(branch.parameter,pass_on,'pass_on');
end
