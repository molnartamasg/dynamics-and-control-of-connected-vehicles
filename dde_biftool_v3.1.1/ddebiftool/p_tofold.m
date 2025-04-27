function fold=p_tofold(funcs,point)
%% convert stst point to fold point
% function fold_point=p_tofold(funcs,point)
% INPUT:
%   funcs problem functions
%	point point near fold point
% OUTPUT:
%	fold_point starting guess for fold point derived from point

% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
% 
%
%%
fold.kind='fold';
fold.parameter=point.parameter;

switch point.kind
    case {'stst','hopf'}
        x=point.x;
    case 'fold'
        error('P_TOFOLD: point is already fold.');
    case 'psol'
        error('P_TOFOLD: periodic psol to fold not supported.');
    otherwise
        if isfield(point,'x')
            x=point.x;
        else
            error('P_TOFOLD: point type %s to fold not supported.',point.kind);
        end
end
fold.x=x;
D=ch_matrix(funcs,x,point.parameter,0);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
fold.v=real(E1(:,i2));
if isfield(point,'stability')
    fold.stability=point.stability;
end
if isfield(point,'nmfm')
    fold.stability=point.nmfm;
end
end
