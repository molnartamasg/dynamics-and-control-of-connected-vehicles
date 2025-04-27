%% Zero-Hopf normal form for zero-Hopf point encountered along fold or Hopf branch
%% Wrapper around nmfm_codimension2_nf
function [nf,nflow,br_ref,indbif]=ZeroHopfNormalform(funcs,branch,inds,varargin)
%% Input
%
% * funcs: problem functions
% * branch: fold or Hopf branch along which zero-Hopf point was encountered
% * inds: array of two successive indices bracing zero-Hopf point
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
% * indbif: is the index in br_ref which is closest to zero-Hopf point
%
% (c) DDE-BIFTOOL v. 3.1.1(72), 30/12/2014
%
%%
type=branch.point(inds(1)).kind;
switch type
    case 'fold'
        detect='hoze';
    case 'hopf'
        detect='zeho';
    otherwise
        error('ZeroHopfNormalform:type',...
            'ZeroHopfNormalform: cannot detect Zero-Hopf along branch of type %s.',...
            type);
end
nmfm_compute=@(f,pt)nmfm_zeho(f,p_tozeho(pt));
[nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin{:});
end