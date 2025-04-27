%% Plot string instability and time to collision indices
clear; %close all; clc;
addpath('data');

%% Create stability charts to plot into
StabilityChartChain;

%% Load results
load('index_maps.mat');
alphat=Ctmap(:,1);  betat=Ctmap(:,2);   Ct=Ctmap(:,3);
alphas=Csmap(:,1);  betas=Csmap(:,2);   Cs=Csmap(:,3);

%% Plot results
% plot ranges
betamin=0;
betamax=1.2;
alphamin=0;
alphamax=1.2;
aspect_ratio=1;

% plot time to collision index
figure(11); hold on; box on;
dotsize = 3.5*5e3;
scatter(betat,alphat,(Ct-0.95*min(Ct)+0.05*max(Ct))*dotsize,'g','filled')
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title('Probability-of-conflict index - C_T');

% plot string instability index
figure(12); hold on; box on;
dotsize = 3.5*700;
scatter(betas,alphas,(Cs-0*min(Cs)+0.05*max(Cs))*dotsize,'r','filled')
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title('String instability index C_S');