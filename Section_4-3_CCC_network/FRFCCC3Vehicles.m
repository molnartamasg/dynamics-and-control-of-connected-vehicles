%% Create frequency response plots for CCC in a 3 vehicle scenario
clear; %close all; clc;

%% Calculate frequency response
% parameters
% HV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;
% CAV
kappa=0.6;
alpha=0.4;
sigma=0.6;
beta1=0.5;
beta2=0.5;  % 0.0, 0.5

% range of frequencies
ommin=0;
ommax=pi;
dom=(ommax-ommin)/1000;
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.5;

% link transfer function
% HV
T12=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
% CAV
D=(-om.^2.*exp(1i*om*sigma)+1i*(alpha+beta1+beta2).*om+alpha*kappa);
T01=(1i*beta1.*om+alpha*kappa)./D;
T02=1i*beta2.*om./D;

% head-to-tail transfer function
G02=T01.*T12+T02;

% get corresponding frequency response
M12=abs(T12);
M01=abs(T01);
M02=abs(T02);
MG=abs(G02);

%% Plot frequency response
figure(1); clf; hold on; box on;
plot(om,M12,om,M01,om,M02,om,MG);
plot([ommin,ommax],[1,1],'k--');
xlabel('\omega [rad/s]');
ylabel('|T(j\omega)|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of CCC with 3 vehicles',10,...
       'HV: \kappa_{\rm h}=',num2str(kappah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
       'CAV: \kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]   ',...
       '\sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       '\beta_1=',num2str(beta1,'%3.2f'),' [1/s]   ',...
       '\beta_2=',num2str(beta2,'%3.2f'),' [1/s]']);
legend('HV: T_{12}','CAV: T_{01}',...
       'CAV: T_{02}','head-to-tail: G_{02}');