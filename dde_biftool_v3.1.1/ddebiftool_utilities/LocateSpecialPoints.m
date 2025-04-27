%% Detect special points along branch of points and compute all normal forms
%% Input
%
% * funcs: problem functions
% * branch: branch of solutions (at the moment supported: hopf, fold, stst)
% * optional name-value pairs:
% * 'max_nferror' (default=|1e-3|) warn if difference between different
% approximation orders is larger than this (useful for finite-difference
% approximation only)
% * all other name-value pairs are passed on as fields for branch fields
%% Output
% * branch: updated branch with inserted special points
% * testfuncs: structure containing fields with branch specific test functions
%
% If an entry in bifpoints is empty detection or computation have failed.
%%
function [branch,testfuncs]=LocateSpecialPoints(funcs,branch,varargin)
%
% (c) DDE-BIFTOOL v. 3.1.1(113), 02/09/2015
%
%% Simple wrapper around branch specific functions
types=struct(...
    'stst',@StstCodimension1,...
    'fold',@FoldCodimension2,...
    'hopf',@HopfCodimension2);
btype=branch.point(1).kind;
routine=types.(btype);
[bifpoints,indices,branch,testfuncs]=routine(funcs,branch,varargin{:});
%% insert bifpoints into refined branch
%#ok<*SFLD>
branch=br_flag(branch);
for i=1:length(indices)
    branch.point(indices(i)).flag=bifpoints{i}.kind;
    branch.point(indices(i)).nmfm=bifpoints{i}.nmfm;
    if isfield(bifpoints{i},'nvec')
        branch.point(indices(i)).nvec=bifpoints{i}.nvec;
    end
end
end