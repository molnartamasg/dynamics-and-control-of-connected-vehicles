%%  DDE-BIFTOOL demo 1 - Neuron: run all scripts of neuron example for testing
% This script runs all demos files in a single run for testing purposes.
% See <html/demo1.html> for the published version of the single scripts
%
% (c) DDE-BIFTOOL v. 3.1.1(125), 29/05/2016
%
%% Description and load path
demo1;
%% Definition of user functions
demo1_funcs;
%% Equilibria
demo1_stst;
%% Hopf bifurcations
demo1_hopf;
%% Normal forms along Hopf bifurcation (broken for r109)
% demo1_normalforms;
%% Periodic orbits
demo1_psol;
%% Homoclinic connection
demo1_hcli;
%% Saddle-node of periodic orbits
demo1_POfold;
%% repeat procedure using convenience functions
demo1_simple;
