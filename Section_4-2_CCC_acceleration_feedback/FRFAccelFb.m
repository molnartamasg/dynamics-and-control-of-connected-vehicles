%% Create frequency response plots for CCC with acceleration feedback
clear; %close all; clc;

%% Calculate frequency response
% parameters
kappa=0.6;
sigma=0.6;
gamma=0.50;	%0.00, 0.25, 0.50, 0.75
alpha=0.4;
beta=0.25;

% range of frequencies
ommin=0;
ommax=6*pi;             %  pi, 6*pi
dom=(ommax-ommin)/600;  % 100,  600
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.2;

% transfer function
T=(alpha*kappa+beta*1i.*om-gamma*om.^2)./...
  ((-om.^2).*exp(1i*om*sigma)+(beta+alpha)*1i.*om+alpha*kappa);
M=abs(T);

%% Plot frequency response
figure(1); clf; hold on; box on;
title([]);
plot([ommin,ommax],[1,1],'k--');
plot([ommin,ommax],[gamma,gamma],'k--');
plot(om,M,'b');
xlabel('\omega [rad/s]');
ylabel('|T(j\omega)|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of the vehicle chain with acceleration feedback',10,...
       '\kappa=',num2str(kappa,'%3.2f'),' [1/s]',',   ',...
       '\sigma=',num2str(sigma,'%3.2f'),' [s]',',   ',...
       '\gamma=',num2str(gamma,'%3.2f'),',   ',...
       '\beta=',num2str(beta,'%3.2f'),' [1/s]',',   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]']);