%% Create analytical stability charts for the vehicle chain

clear; %close all; clc;

%% Parameters
kappa=0.6;
% range of parameters
betamin=-2;
betamax=3;
alphamin=-1;
alphamax=4;
aspect_ratio=1;

%% Plant stability boundaries
betaP1=[betamin,betamax];
alphaP1=[0,0];
alphaP2=[0,alphamax];
betaP2=-alphaP2;

%% Sring stability boundaries
% omega=0
betaS1=[betamin,betamax];
alphaS1=[0,0];
betaS2=[betamin,betamax];
alphaS2=2*(kappa-betaS2);

%% Stability chart
figure(1); clf; hold on; box on;
% plant stability boundaries
plot(betaP1,alphaP1,'r',betaP2,alphaP2,'r');
% string stability boundaries
plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),' [1/s]']);

%% Plot points on existing stability chart
% hold on; plot([-0.5,-0.4,0.25,0.4,0.5,0.5,0.5],[0.4,0.4,0.4,0.4,0.4,0,-0.1],'bx');
% axis([-0.6,1,-0.4,1.2]);