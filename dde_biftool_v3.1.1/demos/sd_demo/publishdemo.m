%% DDE-BIFTOOL  state-dependent delays sd-demo
% This script publishes all demos files in a single run for testing purposes.
% See <html/sd_demo.html> for the published version of the single scripts
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
%% Description and load path
publish('sd_demo.m');
%% Definition of user functions
opts={'maxOutputLines',20};
publish('sd_demo_funcs.m',opts{:});
%% Equilibria
publish('sd_demo_stst.m',opts{:});
%% Hopf bifurcations
publish('sd_demo_hopf.m',opts{:});
%% Periodic orbits
publish('sd_demo_psol.m',opts{:});

