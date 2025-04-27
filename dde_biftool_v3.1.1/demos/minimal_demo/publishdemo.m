%%  DDE-BIFTOOL minimal_demo - Duffing oscillator with delayed feedback
% run all scripts of this for testing
% This script runs all demos files in a single run for testing purposes.
% See <html/minimal_demo.html> for the published version of the single scripts
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(73), 31/12/2014
% </html>
%% Description, load path and define system
publish('minimal_demo');
%% Steady-state bifurcations and periodic orbits
opts={'maxOutputLines',20};
publish('minimal_demo_stst_psol',opts{:});
%% Normal form coefficients for Hopf bifurcations
publish('minimal_demo_extra_nmfm','maxOutputLines',Inf);
%% local bifurcations of periodic orbits
publish('minimal_demo_extra_psol',opts{:});
%% Plot two-parameter bifurcation diagram
publish('minimal_demo_plot_2dbif');
