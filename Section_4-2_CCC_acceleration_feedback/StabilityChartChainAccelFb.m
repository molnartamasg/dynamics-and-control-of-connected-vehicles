%% Create analytical stability charts for CCC with acceleration feedback
clear; %close all; clc;

%% Parameters
kappa=0.6;
sigma=0.6;
gamma=0.25;	% 0, 0.25, 0.5, 0.75
% range of parameters
betamin=-2;         % -2, -0.5
betamax=3;          %  3,  1.5
alphamin=-1;        % -1,    0
alphamax=4;         %  4,    2
aspect_ratio=1;     %  1,    1

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
alphaP2=Om.^2/kappa.*cos(Om*sigma);
betaP2=Om.*sin(Om*sigma)-alphaP2;

% expressions of the om=0 string stability boundaries
betaS1=[betamin,betamax];
alphaS1=[0,0];
betaS2=[betamin,betamax];
alphaS2=2*(kappa*(1-gamma)-betaS2);

% expressions of om>0 string stability boundaries
denab=(2*kappa*sigma-1)*sin(om*sigma)-om*sigma.*cos(om*sigma);
aa=(om*(kappa*sigma-1+gamma^2)+kappa*sin(om*sigma).*(cos(om*sigma)-gamma)...
    -gamma*kappa*om*sigma.*cos(om*sigma))./denab;
bb=((1-gamma^2)*om.^2.*(sin(om*sigma)-om*sigma.*cos(om*sigma)))./denab;
alphaS3=aa-sqrt(aa.^2+bb);
alphaS3(aa.^2+bb<0)=nan;
betaS3=((1-gamma^2)*om+alphaS3*kappa*sigma.*sin(om*sigma))./...
    (sin(om*sigma)+om*sigma.*cos(om*sigma))-alphaS3;
betaS3(betaS3>betamax | betaS3<betamin)=nan; % remove asymptote
alphaS4=aa+sqrt(aa.^2+bb);
alphaS4(aa.^2+bb<0)=nan;
betaS4=((1-gamma^2)*om+alphaS4*kappa*sigma.*sin(om*sigma))./...
    (sin(om*sigma)+om*sigma.*cos(om*sigma))-alphaS4;
betaS4(betaS4>betamax | betaS4<betamin)=nan; % remove asymptote

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
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with acceleration feedback',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \sigma = ',num2str(sigma,'%3.2f'),' [s]',...
       ' [1/s],   \gamma = ',num2str(gamma,'%3.2f')]);