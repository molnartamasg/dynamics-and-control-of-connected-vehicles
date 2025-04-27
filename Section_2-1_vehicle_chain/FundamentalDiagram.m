%% Plot fundamental diagrams
clear; clc; % close all;

% parameters
hst=5;
hgo=55;
vmax=30;
kappa=0.6;  % 1/t_h for IDM; 0.6, 1.2
l=5;
rhost=1/(hst+l);
rhogo=1/(hgo+l);

% plot range
qplotmin=0;
qplotmax=0.75;
rhoplotmin=0;
rhoplotmax=0.12;
drho=rhogo/100;
rho=rhoplotmin:drho:rhoplotmax;

% formulas for IDM
syms v hh k hstop vm
idm=solve(((hstop+v/k)./sqrt(1-(v/vm).^4))==hh,v,'MaxDegree',4);

% range policies
Vlin=@(h)(hst<h & h<hgo).*(h-hst)/(hgo-hst)*vmax + (hgo<=h)*vmax;
Vquad=@(h)(hst<h & h<hgo).*(2*hgo-hst-h).*(h-hst)/(hgo-hst)^2*vmax + (hgo<=h)*vmax;
Vcube=@(h)(hst<h & h<hgo).*(3*hgo-hst-2*h).*(h-hst).^2/(hgo-hst)^3*vmax + (hgo<=h)*vmax;
Vidm=@(hh)eval(subs(idm(2),[k,hstop,vm],[kappa,hst,vmax]));

% fundamental diagram
Qlin=@(rho)(rho>0).*rho.*Vlin(1./rho-l);
Qquad=@(rho)(rho>0).*rho.*Vquad(1./rho-l);
Qcube=@(rho)(rho>0).*rho.*Vcube(1./rho-l);
Qidm=@(rho)(rho>0).*rho.*Vidm(1./rho-l);

% plot fundamental diagrams
figure(1); clf; hold on; box on;
plot([rhost,rhost],[qplotmin,qplotmax],'k--');
plot([rhogo,rhogo],[qplotmin,qplotmax],'k--');
plot([rhoplotmin,rhoplotmax],[rhogo*vmax,rhogo*vmax],'k--');
plot(rho,Qlin(rho),'b');
xlabel('density, \rho (1/m)'); ylabel('flux, Q (1/s)');
axis([rhoplotmin,rhoplotmax,qplotmin,qplotmax]);

figure(2); clf; hold on; box on;
plot([rhost,rhost],[qplotmin,qplotmax],'k--');
plot([rhogo,rhogo],[qplotmin,qplotmax],'k--');
plot([rhoplotmin,rhoplotmax],[rhogo*vmax,rhogo*vmax],'k--');
plot(rho,Qquad(rho),'b');
xlabel('density, \rho (1/m)'); ylabel('flux, Q (1/s)');
axis([rhoplotmin,rhoplotmax,qplotmin,qplotmax]);

figure(3); clf; hold on; box on;
plot([rhost,rhost],[qplotmin,qplotmax],'k--');
plot([rhogo,rhogo],[qplotmin,qplotmax],'k--');
plot([rhoplotmin,rhoplotmax],[rhogo*vmax,rhogo*vmax],'k--');
plot(rho,Qcube(rho),'b');
xlabel('density, \rho (1/m)'); ylabel('flux, Q (1/s)');
axis([rhoplotmin,rhoplotmax,qplotmin,qplotmax]);

figure(4); clf; hold on; box on;
qqplotmin=qplotmin-0.15; qqplotmax=qplotmax-0.15;
plot([rhost,rhost],[qqplotmin,qqplotmax],'k--');
plot([rhogo,rhogo],[qqplotmin,qqplotmax],'k--');
plot([rhoplotmin,rhoplotmax],[0,0],'k--');
plot([rhoplotmin,rhoplotmax],[rhogo*vmax,rhogo*vmax],'k--');
plot(rho,Qidm(rho),'b');
xlabel('density, \rho (1/m)'); ylabel('flux, Q (1/s)');
axis([rhoplotmin,rhoplotmax,qqplotmin,qqplotmax]);