%% Demo - System with cusp and Takens-Bogdanov bifurcations
%
% Demo contributed by Maikel Bosschaert
% 
% from Giannakopoulos, F. and Zapp, A. (2001). Bifurcations in a planar
% system of differential delay equations modeling neural activity. Physica
% D: Nonlinear Phenomena, 159(3):215-232.
%
%% Differential equations
%
% 
% $$\mu\dot{u}_{1}(t)=-u_{1}(t)+q_{11}\alpha(u_{1}(t-T))-q_{12}u_{2}(t-T)+e_{1}$$
%
% $$\mu\dot{u}_{2}(t)=-u_{2}(t)+q_{21}\alpha(u_{1}(t-T))-q_{22}u_{2}(t-T)+e_{2}$$
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(127), 05/09/2016
% </html>
%%
clear % clear variables
close all; % close figures
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');

disp('Cusp demo');

funcs=set_funcs(...
    'sys_rhs', @cusp_rhs,...
    'sys_tau', @()6,...
    'sys_deri', @cusp_deri,...
    'sys_mfderi',@(xx,par,varargin)cusp_mfderi(xx,par,varargin{:}));
%% continue stst in E
Q=0;
E=0;

q11=2.6; q12=Q; q21=1; e1=E; e2=0; T=1;
par = [q11,q12,q21,e1,e2,T];

stst_br = SetupStst(funcs,'x',[0;0],'parameter',par,...
    'contpar',4,'max_step',[4 0.02],'max_bound',[4 1],'min_bound',[4 0],...
    'newheuristics_tests',0);

stst_br.method.continuation.plot = 1;
[stst_br,s,f,r] = br_contn(funcs,stst_br,300);

stst_br.method.bifurcation.plot_testfunctions=1;

stst_br.method.stability.minimal_real_part=-5;
stst_br=br_stabl(funcs,stst_br,0,0);
stst_br=br_bifdet(funcs,stst_br);

%% continue fold point in (E,Q)
FPI=br_getflags(stst_br);
start_ind = FPI(bif2num('fold'),1);

[fold_branch, suc] = SetupFold(funcs, stst_br, start_ind, 'contpar', [2 4], 'dir', 4, 'step', 0.02);

Qmin=-0.1; Qmax=2; Emin=-0.6; Emax=0.6;
fold_branch.parameter.min_bound=[2 Qmin; 4 Emin];
fold_branch.parameter.max_bound=[2 Qmax; 4 Emax];
fold_branch.parameter.max_step=[2 0.02; 4 0.02];

figure;
[fold_branch,s,f,r]=br_contn(funcs,fold_branch,300);
fold_branch = br_rvers(fold_branch);
[fold_branch,s,f,r]=br_contn(funcs,fold_branch,300);

fold_branch.method.bifurcation.plot_testfunctions=1;

fold_branch = br_stabl(funcs,fold_branch,0,0);
fold_branch = br_bifdet(funcs,fold_branch);

