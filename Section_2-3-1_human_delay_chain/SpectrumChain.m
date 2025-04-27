%% Create spectrum plots for the vehicle chain
clear; %close all; clc;

% add ddebiftool folder to path
addpath(genpath('..\dde_biftool_v3.1.1'));

%% Parameters and system definition
% following vehicle
tau=0.6;    % 0.0, 0.6
beta=-0.5;	%  -0.5, 0.5,  0.5
alpha=0.4;	%   0.4, 0.4, -0.1
hst=5;
vmax=30;
kappa=0.6;
amin=7;
amax=3;
% leading vehicle (constant speed for plant stability analysis)
vstar=15;

% range of plot
Remin=-8;   % -1,  -6,  -8
Remax=1;    %  1,   1,   1
Immin=-15;	% -1, -15, -15
Immax=15;	%  1,  15,  15

%% DDE-Biftool for calculating the spectrum
% vector of parameters
par=[vstar,alpha,beta,tau,hst,vmax,kappa,amin,amax];
ind_tau=4;    % index of delay in parameters

% set up right-hand side of system and delays
funcs=set_funcs('sys_rhs',@sys_rhs_chain,'sys_tau',@()ind_tau);

% initial guess for steady state headway, velocity
hguess=30;
vguess=15;

% parameters and initial guess for steady state
stst.kind='stst';
stst.parameter=par;
stst.x=[vguess;hguess];

% correct initial guess for steady state
method=df_mthod(funcs,'stst');
method.stability.minimal_real_part=-8;
% comment this line for alpha=0 and give the correct vguess, hguess
[stst,~]=p_correc(funcs,stst,[],[],method.point);

% compute stability of steady state
stst.stability=p_stabil(funcs,stst,method.stability);

%% Spectrum plot
% plot eigenvalues
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
p_splot(stst);
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(s)');
ylabel('Im(s)');
title(['Spectrum of the vehicle chain',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \tau=',num2str(tau,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);