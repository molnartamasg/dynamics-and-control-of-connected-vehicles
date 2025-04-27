%% Create spectrum plots for the ring configuration
clear; %close all; clc;

% add ddebiftool folder to path
addpath(genpath('..\dde_biftool_v3.1.1'));

%% Parameters and system definition
% number of vehicles
NN=24;      % 24, 8
% vehicle parameters
tau=0.6;    %   0,   0,     0,  0.6, 0.6, 0.6
beta=0.5;	% 0.25, 0.5, 0.25, 0.25, 0.5, 0.8
alpha=0.4;	%  0.4, 0.4, -0.1,  0.4, 0.4, 0.4
hst=5;
vmax=30;
kappa=0.6;
amin=7;
amax=3;
% equilibrium speed
vstar=15;
hstar=hst+vstar/kappa;
% ring length
L=NN*hstar;

% range of plot
Remin=-0.8;	%   -1, -5,  -6, -0.8
Remax=0.2;	%  0.2,  1,   1,  0.2
Immin=-2;	% -0.6, -3, -15,   -2
Immax=2;	%  0.6,  3,  15,    2

%% DDE-Biftool for calculating the spectrum
% vector of parameters
parinit=[L,alpha,beta,tau,hst,vmax,kappa,amin,amax,NN];
ind_tau=4;    % index of delay in parameters

% set up right-hand side of system and delays
funcs=set_funcs('sys_rhs',@sys_rhs_ring,'sys_tau',@()ind_tau);

% initial guess for steady state headway, velocity
vguess=repmat(vstar,NN,1);
hguess=repmat(hstar+1,NN-1,1);

% parameters and initial guess for steady state
stst.kind='stst';
stst.parameter=parinit;
stst.x=[vguess;hguess];

% correct initial guess for steady state
method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=-8; 
[stst,~]=p_correc(funcs,stst,[],[],method.point);

% compute stability of steady state
stst.stability=p_stabil(funcs,stst,method.stability);

%% Spectrum plot
% plot eigenvalues
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
plot(0,0,'b*');
p_splot(stst);
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(s)');
ylabel('Im(s)');
title(['Spectrum of the ring configuration',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \tau=',num2str(tau,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);