%% Plot range policies
clear; clc; % close all;

% parameters
hst=5;
hgo=55;
vmax=30;
kappa=0.6;  % 1/t_h for IDM; 0.6, 1.2

% plot range
vplotmin=0;
vplotmax=40;
hplotmin=0;
hplotmax=70;
dh=0.1;
h=hplotmin:dh:hplotmax;

% plot slope of range policy at h*, between [h*-deltah,h*+deltah]
hstar=35;
deltah=20;

% formulas for IDM
syms v hh k hstop vm
idm=solve(((hstop+v/k)./sqrt(1-(v/vm).^4))==hh,v,'MaxDegree',4);
idmp=diff(idm,hh);
kidm=eval(subs(idmp(2),[k,hstop,vm,hh],[kappa,hst,vmax,hstar]));

% range policies
Vlin=@(h)(hst<h & h<hgo).*(h-hst)/(hgo-hst)*vmax + (hgo<=h)*vmax;
Vquad=@(h)(hst<h & h<hgo).*(2*hgo-hst-h).*(h-hst)/(hgo-hst)^2*vmax + (hgo<=h)*vmax;
Vcube=@(h)(hst<h & h<hgo).*(3*hgo-hst-2*h).*(h-hst).^2/(hgo-hst)^3*vmax + (hgo<=h)*vmax;
Vidm=@(hh)eval(subs(idm(2),[k,hstop,vm],[kappa,hst,vmax]));

% tangents to range policies
Vplin=@(h)(h-hstar)*1/(hgo-hst)*vmax + Vlin(hstar);
Vpquad=@(h)(h-hstar)*(-(hstar-hst)+(2*hgo-hst-hstar))/(hgo-hst)^2*vmax + Vquad(hstar);
Vpcube=@(h)(h-hstar)*(-2*(hstar-hst)^2+(3*hgo-hst-2*hstar)*2*(hstar-hst))/(hgo-hst)^3*vmax + Vcube(hstar);
Vpidm=@(h)(h-hstar)*kidm+Vidm(hstar);

% plot range policies
figure(1); clf; hold on; box on;
plot([hst,hst],[vplotmin,vplotmax],'k--');
plot([hgo,hgo],[vplotmin,vplotmax],'k--');
plot([hplotmin,hplotmax],[vmax,vmax],'k--');
plot(h,Vlin(h),'b');
plot([hstar-deltah,hstar+deltah],Vplin([hstar-deltah,hstar+deltah]),'m');
xlabel('headway, h (m)'); ylabel('desired speed, V (m/s)');
axis([hplotmin,hplotmax,vplotmin,vplotmax]);

figure(2); clf; hold on; box on;
plot([hst,hst],[vplotmin,vplotmax],'k--');
plot([hgo,hgo],[vplotmin,vplotmax],'k--');
plot([hplotmin,hplotmax],[vmax,vmax],'k--');
plot(h,Vquad(h),'b');
plot([hstar-deltah,hstar+deltah],Vpquad([hstar-deltah,hstar+deltah]),'m');
xlabel('headway, h (m)'); ylabel('desired speed, V (m/s)');
axis([hplotmin,hplotmax,vplotmin,vplotmax]);

figure(3); clf; hold on; box on;
plot([hst,hst],[vplotmin,vplotmax],'k--');
plot([hgo,hgo],[vplotmin,vplotmax],'k--');
plot([hplotmin,hplotmax],[vmax,vmax],'k--');
plot(h,Vcube(h),'b');
plot([hstar-deltah,hstar+deltah],Vpcube([hstar-deltah,hstar+deltah]),'m');
xlabel('headway, h (m)'); ylabel('desired speed, V (m/s)');
axis([hplotmin,hplotmax,vplotmin,vplotmax]);

figure(4); clf; hold on; box on;
vvplotmin=vplotmin-5; vvplotmax=vplotmax-5;
plot([hst,hst],[vvplotmin,vvplotmax],'k--');
plot([hgo,hgo],[vvplotmin,vvplotmax],'k--');
plot([hplotmin,hplotmax],[0,0],'k--');
plot([hplotmin,hplotmax],[vmax,vmax],'k--');
plot(h,Vidm(h),'b');
plot([hstar-deltah,hstar+deltah],Vpidm([hstar-deltah,hstar+deltah]),'m');
xlabel('headway, h (m)'); ylabel('desired speed, V (m/s)');
axis([hplotmin,hplotmax,vvplotmin,vvplotmax]);