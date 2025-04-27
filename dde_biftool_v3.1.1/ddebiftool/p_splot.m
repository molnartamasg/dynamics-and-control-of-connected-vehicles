%% Plot spectrum of point with stability information
%%
function p_splot(point,varargin)
% function p_splot(point)
% INPUT:
%	point point whose stability needs plotting
%  'plotaxis' (optional, default gca): axis on which to plot
%
% (c) DDE-BIFTOOL v. 3.1.1(124), 05/02/2016
%
%%
default={'plotaxis',gca};
options=dde_set_options(default,varargin,'pass_on');
switch point.kind
    case {'stst','fold','hopf'}
        if isfield(point.stability,'l0')
            root_plt(point.stability.l0,point.stability.l1,point.stability.n1,'plotaxis',options.plotaxis);
        end
        xlabel('\Re(\lambda)');
        ylabel('\Im(\lambda)');
    case 'psol'
        if isfield(point.stability,'mu')
            mult_plt(point.stability.mu,'plotaxis',options.plotaxis);
        end
        xlabel('\Re(\mu)');
        ylabel('\Im(\mu)');
    otherwise
        err=point.kind;
        error('P_SPLOT: point kind %s not recognized.',err);
end
end
