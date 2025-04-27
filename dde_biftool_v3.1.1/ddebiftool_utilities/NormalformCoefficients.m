%% Normal form coefficient a in fold point along fold or Hopf branch
%% Input
%
% * funcs: problem functions
% * branch: fold branch or array of points of kind 'fold' or 'hopf'
% * optional name-value pair: 
% * 'use_nullvectors' (boolean, default false, for hopf only) use nullvectors from
% previous point along branch to compute new nullvectors by bordering
% (otherwise, svd is used)
%
%% Output
%
% * L1: array of normal form coefficients along branch
% * L1low: equals b if funcs.sys_mfderi_provided is true, otherwise (for
% finite-difference approximation) a lower-order estimate of L1
% coefficients. Use L1-L1low to estimate the accuracy of L1. If this is
% large, finite-differences are problematic. If this is small it may(!) be
% ok.
%%
function [L1,L1low]=NormalformCoefficients(funcs,branch,varargin)
%%
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
%%
default={'use_nullvectors',false,'print',0};
options=dde_set_options(default,varargin,'pass_on');
if isfield(branch,'point')
    pt=branch.point;
else
    pt=branch;
end
switch pt(1).kind
    case 'hopf'
        nffunc=@nmfm_hopf;
        cfield='L1';
    case 'fold'
        nffunc=@(fun,point,dummy)nmfm_fold(fun,point);
        cfield='b';
end
npt=length(pt);
nullpoint={};
L1=NaN(1,npt);
L1low=L1;
for i=1:npt
    currpt=pt(i);
    nmfm_compute=@(fun,point)nffunc(fun,point,nullpoint{:});
    [nf,nflow]=nmfm_wrap_findiff(funcs,currpt,nmfm_compute);
    L1(i)=nf.nmfm.(cfield);
    L1low(i)=nflow.nmfm.(cfield);
    if options.print>0
        fprintf('Normalform Coefficients: pt %d of %d on %s branch: %s=%g\n',...
            i,npt,pt(1).kind,cfield,L1(i));
    end
    if options.use_nullvectors
        nullpoint={nf};
    end
end
end
