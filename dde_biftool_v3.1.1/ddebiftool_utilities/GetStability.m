%% Number of unstable eigenvalues/Floquet multipliers along branch
%%
function [nunst,dom,triv_defect,points]=GetStability(branch,varargin)
%% Inputs
%
% * |branch|: branch along which number of unstable eigenvalues/Floquet
% multipliers is determined
%
%  Important name-value pair inputs
% 
% * |'exclude_trivial'| (logical): exclude trivial eigenvalues (what they are
%                            depends on the type of points)
% * |'locate_trivial'| (function): e.g. @(p)[1;-1] for removing 1 and -1
%                            for period doubling orbits (overwrites standard location)
% * |'funcs'|:   set of functions for computing branch, needed if stability 
% information is not yet present and needs to be calculated
% * |'method'|: which method is used for comutation (if not present, a
% standard method is generated with |df_mthod|)
% * |'recompute'| (default |false|):  force recomputation of stability even
% if it is present.
%
%% Outputs
% * |nunst| (vector of integers): number of unstable eigenvalues
% * |dom| (vector of real/complex): dominant eigenvalue (closest unstable to
%                       bifurcation if exists, otherwise, closest stable)
% * |triv_defect| (vector of real): distance of eigenvalue approximating
%                               trivial eigenvalue from its true value
% * |points| (struct array): array of points (same as in input branch but
% with stability info if not present on input)
%
% (c) DDE-BIFTOOL v. 3.1.1(73), 31/12/2014
%
%% process options
defaults={'exclude_trivial',false,'locate_trivial',[],'funcs',[],'critical',false,...
    'points',[],'method',[],'recompute',false,'stabilityfield','l0'};
options=dde_set_options(defaults,varargin);
%% check if stability has been computed already & compute if necessary
if isfield(branch,'point')
    points=branch.point;
    if ~isempty(options.points) % select only some points if desired
        points=branch.point(options.points);
    end
else
    points=branch;
end
if isfield(branch,'method') && isfield(branch.method,'stability')
    mth=branch.method;
elseif ~isempty(options.method)
    mth=options.method;
else
    mth=df_mthod(options.funcs,points(1).kind);
end
for i=1:length(points)
    if options.recompute || ~isfield(points(i),'stability') ||...
            isempty(points(i).stability)
        stab=p_stabil(options.funcs,points(i),mth.stability);
        points(i).stability=stab;
    end
end
%% Count nmber of unstable eigenvalues
% (removing trivial ones as indicated by |locate_trivial| input or standard
% theory)
np=length(points);
nunst=NaN(np,1);
dom=nunst;
triv_defect=zeros(np,1);
switch points(1).kind
    case 'psol'
        getev=@(p)p.stability.mu;
        stab=@(x)log(abs(x)+eps);
        triv=@(p)1;
    case 'stst'
        getev=@(p)p.stability.(options.stabilityfield);
        stab=@(x)real(x);
        triv=@(p)[];
    case 'fold'
        getev=@(p)p.stability.(options.stabilityfield);
        stab=@(x)real(x);
        triv=@(p)0;
    case 'hopf'
        getev=@(p)p.stability.(options.stabilityfield);
        stab=@(x)real(x);
        triv=@(p)[1i*p.omega;-1i*p.omega];
    otherwise
        error('GetStability:type not implemented');
end
if ~isempty(options.locate_trivial)
    % overwrite location of trivial eigenvalues
    triv=options.locate_trivial;
end
for i=1:np
    p=points(i);
    excl=triv(p);
    ev=getev(p);
    if options.exclude_trivial
        for k=1:length(excl)
            [dist,ind]=min(abs(ev-excl(k)));
            ev(ind)=[];
            if ~isempty(dist)
                triv_defect(i)=max(dist,triv_defect(i));
            end
        end
    end
    unstsel=stab(ev)>=0;
    nunstloc=sum(unstsel);
    if isempty(nunstloc)
        nunst(i)=0;
    else
        nunst(i)=nunstloc;
    end
    if nunst(i)>0 && ~options.critical
        [dum,ind]=min(stab(ev(unstsel))); %#ok<ASGLU>
        dom(i)=ev(ind);
        if imag(ev(ind))<0
            dom(i)=conj(dom(i));
        end
    else
        [dum,ind]=min(abs(stab(ev)));
        if ~isempty(dum)
            dom(i)=ev(ind);
            if imag(ev(ind))<0
                dom(i)=conj(dom(i));
            end
        end
    end
end
end
