%% Fold normal form for fold point encountered along stst branch
%% Wrapper around SetupStstBifurcation and nmfm_{fold,hopf}
function [nf,nflow,br_ref,indbif]=Codim1Normalform(funcs,branch,inds,type,varargin)
%% Input
%
% * funcs: problem functions
% * branch: stst branch along which bifurcation point was encountered
% * inds: array of two successive indices bracing bifurcation point
%
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
% * br_ref: refinement using extended system, including sinlge new point
% along branch
% * indbif: is the index in br_ref which is closest to fold point
%
% (c) DDE-BIFTOOL v. 3.1.1(115), 02/09/2015
%
%%
default={'print',0};
options=dde_set_options(default,varargin,'pass_on'); 
[bifbranch,suc]=SetupStstBifurcation(funcs,branch,inds(1),type,varargin{:});
if ~suc
    if options.print>0
        fprintf('Comdim1Normalform: nonlinear system not successfully solved\n');
    end
    nf=[];
    nflow=[];
    br_ref=branch;
    indbif=[];
    return
elseif options.print>0
    fprintf('Codim1Normalform: %s point found\n',type);
end
bifpoint=bifbranch.point(1);
bifstst=p_tostst(funcs,bifpoint);
br_ref=branch;
if isfield(br_ref.point(1),'stability')
    bifstst.stability=p_stabil(funcs,bifstst,branch.method.stability);
end
bifstst=nmfm_trim_point(bifstst,br_ref.point(1));
br_ref.point=[br_ref.point(1:inds(1)),bifstst,br_ref.point(inds(2):end)];
indbif=inds(2);
nmfm_compute=str2func(['nmfm_',type]);
[nf,nflow]=nmfm_wrap_findiff(funcs,bifpoint,nmfm_compute);
end