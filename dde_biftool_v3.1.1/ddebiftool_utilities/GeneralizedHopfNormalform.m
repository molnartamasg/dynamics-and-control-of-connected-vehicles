%% generalized Hopf normal form for generalized Hopf point encountered along Hopf branch
%% Wrapper around nmfm_codimension2_nf
%% Input
%
% * funcs: problem functions
% * branch: Hopf branch along which generalized Hopf point was encountered
% * inds: array of two successive point indices bracing generalized Hopf point
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
% (c) DDE-BIFTOOL v. 3.1.1(72), 30/12/2014
%
%%
function [nf,nflow,br_ref,indbif]=GeneralizedHopfNormalform(funcs,branch,inds,varargin)
if ~strcmp(branch.point(inds(1)).kind,'hopf')
    error('GeneralizedHopfNormalform:type',...
            'GeneralizedHopfNormalform: cannot detect generalized Hopf along branch of type %s.',...
            type);
end
detect=@(p)genhdetect(funcs,p);
nmfm_compute=@(f,pt)nmfm_genh(f,p_togenh(pt));
[nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin{:});
end
function L1=genhdetect(funcs,p)
newpoint=nmfm_hopf(funcs,p);
L1=newpoint.nmfm.L1;
end