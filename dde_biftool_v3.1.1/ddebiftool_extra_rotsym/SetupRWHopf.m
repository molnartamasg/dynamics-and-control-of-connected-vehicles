%% SetupRWHopf - Initialize continuation of Hopf bifurcations of relative equilibria
%%
function [hbranch,suc]=SetupRWHopf(funcs,branch,ind,varargin)
%% 
% simple wrapper to have a sensible name (see demo rotsym_demo how to use this function).
%
% (c) DDE-BIFTOOL v. 3.1.1(29), 15/04/2014
%
[hbranch,suc]=SetupHopf(funcs,branch,ind,varargin{:});
end
