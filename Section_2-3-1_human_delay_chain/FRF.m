%% Create frequency response plots for the vehicle chain with time delay

clear; %close all; clc;

%% Calculate frequency response
% parameters
tau=0.6;	% 0.0, 0.6
kappa=0.6;
alpha=0.4;
beta=0.75143;	% 0.25, 0.40, 0.50, 0.75143, 0.80

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
  ((-om.^2).*exp(1i*om*tau)+(beta+alpha)*1i.*om+alpha*kappa);
M=abs(T);
P=om.^2-2*alpha*kappa*cos(om*tau)-2*(alpha+beta)*om.*sin(om*tau)+alpha*(alpha+2*beta);

%% Plot frequency response
figure(1); clf; hold on; box on;
title([]);
plot(om,M);
plot([ommin,ommax],[1,1],'k--');
xlabel('\omega [rad/s]');
ylabel('|T(j\omega)|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of the vehicle chain with time delay',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \tau=',num2str(tau,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);
   
figure(2); clf; hold on; box on;
title([]);
plot(om,P);
plot([ommin,ommax],[0,0],'k--');
xlabel('\omega [rad/s]');
ylabel('P(\omega)');
axis([ommin,ommax,Pmin,Pmax]);
title(['Frequency response of the vehicle chain with time delay',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \tau=',num2str(tau,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);