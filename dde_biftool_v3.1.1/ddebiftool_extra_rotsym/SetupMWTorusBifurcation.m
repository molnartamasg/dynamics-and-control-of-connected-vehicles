%% SetupMWTorusBifurcation - Initialize continuation of torus bifurcations of relative periodic orbits
%%
function [trfuncs,trbranch,suc]=SetupMWTorusBifurcation(funcs,branch,ind,varargin)
%%
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two (see demo rotsym_demo how to do this).
%
% (c) DDE-BIFTOOL v. 3.1.1(29), 15/04/2014
%
[trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,'nremove',2,varargin{:});
end
