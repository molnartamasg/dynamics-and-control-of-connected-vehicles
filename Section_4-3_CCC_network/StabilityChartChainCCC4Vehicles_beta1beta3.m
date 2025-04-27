%% Create stability charts for CCC in a 4 vehicle scenario
clear; %close all; clc;

%% Parameters
% HV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;
% CAV
kappa=0.6;
alpha=0.4;
sigma=0.6;
beta2=0.5; % 0.0, 0.5, 1.0

% range of parameters
beta1min=-0.5;  % -0.5,   0
beta1max=1.5;   %  1.5, 1.2
beta3min=-0.5;  % -0.5,   0
beta3max=1.5;   %  1.5, 1.2
beta1st=200;
beta3st=200;
aspect_ratio=1;

% range of frequencies
ommin=0;
ommax=2*pi;
dom=2*pi/200;

%% Frequency response
% grid in plane of parameters
beta1_v=linspace(beta1min,beta1max,beta1st+1);	% vector of parameter values
beta3_v=linspace(beta3min,beta3max,beta3st+1);	% vector of parameter values
om_v=ommin:dom:ommax;                           % vector of frequency values
[beta1,beta3,om]=meshgrid(beta1_v,beta3_v,om_v);% matrix of parameter and frequency

% link transfer function
% HV
T12=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
T23=T12;
% CAV
D=(-om.^2.*exp(1i*om*sigma)+1i*(alpha+beta1+beta2+beta3).*om+alpha*kappa);
T01=(1i*beta1.*om+alpha*kappa)./D;
T02=1i*beta2.*om./D;
T03=1i*beta3.*om./D;

% head-to-tail transfer function
G03=T01.*T12.*T23+T02.*T23+T03;

% matrix of maximum head-to-tail transfer function magnitudes
Mmax = max(abs(G03(:,:,0<om_v)),[],3);

% plant stability boundaries
Om1=fsolve(@(Om)alpha-Om^2*cos(Om*sigma)/kappa,0.5);
Om2=fsolve(@(Om)alpha-Om^2*cos(Om*sigma)/kappa,2.5);
sumbetaI=Om1*sin(Om1*sigma)-alpha;
sumbetaII=Om2*sin(Om2*sigma)-alpha;

%% Stability chart
figure(1); clf; hold on; box on;
contourf(beta1(:,:,1),beta3(:,:,1),-Mmax,-[1 1],'b','Linewidth',1.5);
colormap gray
plot([beta1min,beta1max],sumbetaI-beta2-[beta1min,beta1max],'r');
plot([beta1min,beta1max],sumbetaII-beta2-[beta1min,beta1max],'r');
axis([beta1min beta1max beta3min beta3max]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_1 [1/s]');
ylabel('\beta_3 [1/s]');
title(['Stability chart of CCC with 4 vehicles',10,...
       'HV: \kappa_h=',num2str(kappah,'%3.2f'),' [1/s]   ',...
       '\alpha_h=',num2str(alphah,'%3.2f'),' [1/s]   ',...
       '\beta_h=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
       'CAV: \kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]   ',...
       '\beta_{2}=',num2str(beta2,'%3.2f'),' [1/s]   ',...
       '\sigma=',num2str(sigma,'%3.2f'),' [s]']);
   
%% Plot points on existing stability chart
% hold on; plot([0.5,0.5],[0,0.5],'bx');