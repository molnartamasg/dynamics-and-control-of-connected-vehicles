%% Normal form for Cusp point encountered along fold branch
%% Wrapper around nmfm_codimension2_nf
%% Input
%
% * funcs: problem functions
% * branch: fold branch along which cusp point was encountered
% * inds: array of two successive point indices bracing cusp point
%
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
% * br_ref: bisection refinements are performed along the branch before
% normal form calculation, br_ref is branch with additional points between
% indices inds
% * indbif: is the index in br_ref which is closest to generalized Hopf point
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
%%
function [nf,nflow,br_ref,indbif]=CuspNormalform(funcs,branch,inds,varargin)
if ~strcmp(branch.point(inds(1)).kind,'fold')
    error('CuspNormalform:type',...
            'CuspNormalform: cannot detect cusp along branch of type %s.',...
            type);
end
detect=@(p)cuspdetect(funcs,p);
nmfm_compute=@(f,pt)nmfm_cusp(f,p_tocusp(pt));
[nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin{:});
end
function L1=cuspdetect(funcs,p)
newpoint=nmfm_fold(funcs,p);
L1=newpoint.nmfm.b;
end