%%  plot Floquet multipliers of linear problem (periodic orbits)
%%
function mult_plt(mu,varargin)
%% function mult_plt(mu)
% INPUT:
%	mu approximated multipliers
%   'plotaxis' (named optional, default: gca) axis on which to plot
%
% Different colors for unstable (red), stable (green) and critical(black)
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
%% process options
default={'plotaxis',gca};
options=dde_set_options(default,varargin,'pass_on');

theta=0:0.25:360;

rmu=mu(abs(mu)>1);
gmu=mu(abs(mu)<1);
kmu=mu(abs(mu)==1);

plot(options.plotaxis,...
    cosd(theta),sind(theta),'b-',...
    real(rmu),imag(rmu),'rx',...
    real(gmu),imag(gmu),'gx',...
    real(kmu),imag(kmu),'kx');
end
