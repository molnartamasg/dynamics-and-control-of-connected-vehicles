%% Create frequency response plots for ACC with time delay
clear; %close all; clc;

%% Calculate frequency response
% number of vehicles
N=12;               % total
% penetration of automation
p=(0:N).'/N;        % AV
% plotting frequency response only for selected penetrations 
psel=(0:2:N).'/N;

% type of vehicle
vehtype='car';
% vehtype='truck';

% vehicle parameters
if strcmp(vehtype,'car')	% personal vehicle
    k=0.45;     % air drag
    m=1500;     % mass
else	% truck
    k=3.85;     % air drag
    m=29500;	% mass
end

% leader speed
vstar=20;
% linear coefficient
d=2*k/m*vstar;

% parameters
% HV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;    % 0.9, 1.0
% CAV
kappa=0.6;
sigma=0.6;
alpha=0.4;
beta=0.5;
eta=0.05;

% range of frequencies
ommin=0;
ommax=pi;
dom=(ommax-ommin)/1000;
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.5;
Mmaxmin=0.7;
Mmaxmax=1.3;

% link transfer function
% HV
Th=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
% CAV
T=(-beta*om.^2+1i*alpha*kappa*om+eta*kappa)./...
  ((-1i*om.^3-d*om.^2).*exp(1i*om*sigma)-(alpha+beta)*om.^2+1i*(alpha*kappa+eta)*om+eta*kappa);

% normalized head-to-tail transfer function
G=T.^p.*Th.^(1-p);

% get corresponding frequency response
Mh=abs(Th);
M=abs(T);
MG=abs(G);

% maximum of frequency response
% MGmax=max(MG,[],2);
[~,idx]=max(Mh,[],2);
MGmax=MG(:,idx);

%% Plot frequency response
figure(1); clf; hold on; box on;
plot(om,Mh,om,M,om,MG(ismember(p,psel),:));
plot([ommin,ommax],[1,1],'k--');
plot([om(idx),om(idx)],[0,Mh(idx)],'k--');
xlabel('\omega [rad/s]');
ylabel('|T(j\omega)|');
axis([ommin,ommax,Mmin,Mmax]);
title(['Frequency response of ',num2str(N),' vehicles',10,...
       'HV: \kappa_{\rm h}=',num2str(kappah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
       'CAV: \kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]   ',...
       '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
       '\eta=',num2str(eta,'%3.2f'),' [1/s^2]   ',...
       '\sigma=',num2str(sigma,'%3.2f'),' [s]']);
legend(horzcat({'HV: T_{\rm h}','CAV: T'},...
    cellstr(num2str(psel*100', 'H2T: G - %.0f%% AV')).'));

% plot maximum of frequency response vs AV penetration
figure(2); clf; hold on; box on;
plot(p*100,MGmax);
plot([0,100],[1,1],'k--');
xlabel('AV penetration [%]');
ylabel('|G(j\omega*)|');
axis([0,100,Mmaxmin,Mmaxmax]);
pbaspect([1,1/3,1]);
title(['Magnification at \omega*=',num2str(om(idx)),' [rad/s]']);