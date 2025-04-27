%% Create analytical stability charts for the vehicle chain with lag and delay

clear; %close all; clc;

%% Parameters
kappa=0.6;
xi=0.6;     % 0.5, 0.6, 0.7, 0.60, 0.48, 0.36, 0.24, 0.12, 0.00    % lag
tau=0.0;	% 0.0, 0.0, 0.0, 0.00, 0.12, 0.24, 0.36, 0.48, 0.60    % delay
% range of parameters
betamin=-0.6;	% -2,  -2,   0, -0.6
betamax=1;      %  3,   4, 1.2,    1
alphamin=-0.4;	% -1,   0,   0, -0.4
alphamax=1.2;	%  4,   4, 1.2
aspect_ratio=1;	%  1, 1.5,   1,  1.2

% range of frequencies for plant stability boundaries
Ommin=0;
Ommax=2*pi;
dOm=2*pi/200;
Om=Ommin:dOm:Ommax;
% range of frequencies for string stability boundaries
ommin=0;
ommax=2*pi;
dom=2*pi/10000;
om=ommin:dom:ommax;

%% Stability boundaries
% expression of the s=0 plant stability boundary
betaP1=[betamin,betamax];
alphaP1=[0,0];

% expression of the s=jOm plant stability boundary
alphaP2=1/kappa*(-Om.^3*xi.*sin(Om*tau)+Om.^2.*cos(Om*tau));
betaP2=Om.^2*xi.*cos(Om*tau)+Om.*sin(Om*tau)-alphaP2;

% expressions of the om=0 string stability boundaries
betaS1=[betamin,betamax];
alphaS1=[0,0];
betaS2=[betamin,betamax];
alphaS2=2*(kappa-betaS2);

% expressions of om>0 string stability boundaries
dencd=(1-om.^2*xi*tau).*sin(om*tau)+om*(tau+2*xi).*cos(om*tau);
cc=om.*(1+2*om.^2*xi^2)./dencd;
dd=kappa*((tau+xi)*sin(om*tau)+om*xi*tau.*cos(om*tau))./dencd-1;
aa=((dd+1-kappa*xi).*om.*sin(om*tau)+((dd+1).*om.^2*xi+kappa).*cos(om*tau)-cc)./(2*dd+1);
bb=(2*cc.*om.*(sin(om*tau)+om*xi.*cos(om*tau))-om.^2-om.^4*xi^2)./(2*dd+1);
alphaS3=aa-sqrt(aa.^2+bb);
alphaS3(aa.^2+bb<0)=nan;
betaS3=cc+dd.*alphaS3;
betaS3(betaS3>betamax | betaS3<betamin)=nan; % remove asymptote
alphaS4=aa+sqrt(aa.^2+bb);
alphaS4(aa.^2+bb<0)=nan;
betaS4=cc+dd.*alphaS4;
betaS4(betaS4>betamax | betaS4<betamin)=nan; % remove asymptote

% parabola corresponding to the om>0 string stability boundaries
if tau==0
    betaS5=linspace(betamin,betamax,100);
    alphaS5=(betaS5-1/2/xi).^2/(1/xi-2*kappa);
end

%% Stability chart
figure(1); clf; hold on; box on;
% plot s=0 plant stability boundary
plot(betaP1,alphaP1,'r');
% plot s=jOm plant stability boundary
plot(betaP2,alphaP2,'r');
% plot om=0 string stability boundaries
plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
% plot om>0 string stability boundaries
plot(betaS3,alphaS3,'b',betaS4,alphaS4,'b');
% plot parabola corresponding to om>0 string stability boundaries
if tau==0
    plot(betaS5,alphaS5,'b--');
end
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with lag and time delay',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \tau = ',num2str(tau,'%3.2f'),' [s]',...
       ' [1/s],   \xi = ',num2str(xi,'%3.2f'),' [s]']);