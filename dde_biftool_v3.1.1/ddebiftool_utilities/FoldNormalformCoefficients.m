%% Normal form coefficient a in fold point along fold branch 
% (wrapper for NormalformCoefficients)
%% Input
%
% * funcs: problem functions
% * branch: fold branch or array of points of kind 'fold'
% * optional name-value pair: 
% * 'sys_mfderi': update sys_mfderi field of funcs (useful for
% finite-difference approximation only)
%
%% Output
%
% * a: array of normal form coefficients along branch
% * a_low: equals b if funcs.sys_mfderi_provided is true, otherwise (for
% finite-difference approximation) a low-er-order estimate of b
% coefficients. Use L1-L1low to estimate the accuracy of b. If this is
% large, finite-differences are problematic. If this is small it may(!) be
% ok.
%%
function [a,a_low]=FoldNormalformCoefficients(funcs,branch,varargin)
%%
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
%%
[a,a_low]=NormalformCoefficients(funcs,branch,varargin{:});
end
