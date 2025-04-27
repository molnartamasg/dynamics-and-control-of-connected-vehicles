%% Create spectrum plots for the ring configuration with no delay
clear; %close all; clc;

%% Parameters and spectrum calculation
% number of vehicles
NN=8;       % 24, 8
% vehicle parameters
beta=0.5;	% 0.25, 0.5, 0.25
alpha=0.4;	%  0.4, 0.4, -0.1
kappa=0.6;

% range of plot
Remin=-1;
Remax=0.2;
Immin=-0.6;
Immax=0.6;

% calculation of spectrum
kk=0:NN-1;
b=alpha+beta*(1-exp(1i*2*kk*pi/NN));
c=alpha*kappa*(1-exp(1i*2*kk*pi/NN));
lambda=[-b/2-sqrt(b.^2-4*c)/2,-b/2+sqrt(b.^2-4*c)/2];

%% Spectrum plot
figure(2); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
plot(real(lambda),imag(lambda),'b*');
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title(['Spectrum of the ring configuration with no delay',10,...
       '\kappa=',num2str(kappa,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);