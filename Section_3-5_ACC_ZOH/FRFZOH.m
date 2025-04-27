%% Create frequency response plots for the vehicle chain with ZOH
clear; %close all; clc;

%% Calculate frequency response
% parameters
kappa=0.6;
deltat=0.4;
alpha=0.4;
beta=0.8;	% 0.25, 0.50, 0.80

% range of frequencies
ommin=0;
ommax=4*pi/deltat;      % 4*pi/deltat, pi
dom=(ommax-ommin)/1600;	% /1600, /100
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.2;

% transfer function
z=exp(1i*om*deltat);
T=((deltat^2/2*alpha*kappa+deltat*beta)*z+deltat^2/2*alpha*kappa-deltat*beta)./...
    (z.^3-2*z.^2+(deltat^2/2*alpha*kappa+deltat*(alpha+beta)+1)*z+deltat^2/2*alpha*kappa-deltat*(alpha+beta));
M=abs(T);

% transfer function of the continuous-time system with same average delay
tau=1.5*deltat;
Tconti=(alpha*kappa+beta*1i.*om)./...
  ((-om.^2).*exp(1i*om*tau)+(beta+alpha)*1i.*om+alpha*kappa);
Mconti=abs(Tconti);

%% Plot frequency response
figure(1); clf; hold on; box on;
plot([ommin,ommax],[1,1],'k--');
plot([2*pi/deltat,2*pi/deltat],[Mmin,Mmax],'k--');
plot(om,Mconti,'b--');
plot(om,M,'b');
xlabel('\omega [rad/s]');
ylabel('|T(e^{j\omega\Deltat})|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of the vehicle chain with ZOH',10,...
       '\kappa=',num2str(kappa,'%3.2f'),' [1/s]',...
       '   \Deltat=',num2str(deltat,'%3.2f'),' [s]',...
       '   \beta=',num2str(beta,'%3.2f'),' [1/s]',...
       '   \alpha=',num2str(alpha,'%3.2f'),' [1/s]']);