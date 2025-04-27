%% Create frequency response plots for CCC with time delay
clear; %close all; clc;

%% Calculate frequency response
% number of vehicles
N=12;               % total
% penetration of automation
p=[0,2,4,6].'/N;	% CAV
% plotting frequency response only for selected penetrations 
psel=p;

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
betan=0.5;  % 0.0, 0.1, 0.2, 0.3, 0.4, 0.5

% range of frequencies
ommin=0;
ommax=pi;
dom=(ommax-ommin)/1000;
om=ommin:dom:ommax;

% range of magnitude on frequency response plot
Mmin=0;
Mmax=1.5;
% Mmaxmin=0.7;
% Mmaxmax=1.3;

% link transfer function
% HV
Th=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
% CAV
D=-om.^2.*exp(1i*om*sigma)+1i*(alpha+beta1+betan).*om+alpha*kappa;
T=(1i*beta1.*om+alpha*kappa)./D;
Tn=1i*betan.*om./D;

% normalized head-to-tail transfer function
G=(T.*Th.^(1./p-1)+Tn).^p;
G(p==0,:)=Th;
% G(p==1,:)=T;

% get corresponding frequency response
Mh=abs(Th);
M=abs(T);
Mn=abs(Tn);
MG=abs(G);

% maximum of frequency response
% MGmax=max(MG,[],2);
[~,idx]=max(Mh,[],2);
MGmax=MG(:,idx);

%% Plot frequency response
figure(1); clf; hold on; box on;
plot(om,Mh,om,M,om,Mn,om,MG(ismember(p,psel),:));
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
       '\beta_1=',num2str(beta1,'%3.2f'),' [1/s]   ',...
       '\beta_n=',num2str(betan,'%3.2f'),' [1/s]   ',...
       '\sigma=',num2str(sigma,'%3.2f'),' [s]']);
legend(horzcat({'HV: T_{\rm h}','CAV: T','CAV: T_n'},...
    cellstr(num2str(psel*100', 'H2T: G - %.0f%% CAV')).'));

% % plot maximum of frequency response vs CAV penetration
% figure(2); clf; hold on; box on;
% plot(p*100,MGmax);
% plot([0,100],[1,1],'k--');
% xlabel('CAV penetration [%]');
% ylabel('|G(j\omega*)|');
% axis([0,100,Mmaxmin,Mmaxmax]);
% pbaspect([1,1/3,1]);
% title(['Magnification at \omega*=',num2str(om(idx)),' [rad/s]']);