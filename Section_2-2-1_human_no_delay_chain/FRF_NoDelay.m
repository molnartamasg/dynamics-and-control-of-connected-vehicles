%% Create frequency response plots for the vehicle chain with time delay

clear; %close all; clc;

%% Calculate frequency response
% parameters
kappa=0.6;
alpha=0.4;
beta=0.25;	% 0.25, 0.40, 0.50

% range of frequencies
ommin=0;
ommax=pi;
dom=(ommax-ommin)/100;
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.2;
Pmin=-0.5;
Pmax=2;

% transfer function
T=(alpha*kappa+beta*1i.*om)./...
  (-om.^2+(beta+alpha)*1i.*om+alpha*kappa);
M=abs(T);

%% Plot frequency response
figure(1); clf; hold on; box on;
plot(om,M);
plot([ommin,ommax],[1,1],'k--');
xlabel('\omega [rad/s]');
ylabel('|T(j\omega)|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of the vehicle chain',10,...
       '\kappa=',num2str(kappa,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);