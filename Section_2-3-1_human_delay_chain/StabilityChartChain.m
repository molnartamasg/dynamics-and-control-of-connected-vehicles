%% Create analytical stability charts for the vehicle chain with time delay
clear; %close all; clc;

%% Parameters
kappa=0.6;	% 0.4, 0.6, 0.8, 1
tau=0.6;	% 0.0, 0.5, 0.6, 0.7
% range of parameters
betamin=-2;         % -2,  -2,   0
betamax=3;          %  3,   4, 1.2
alphamin=-1;        % -1,   0,   0
alphamax=4;         %  4,   4, 1.2
aspect_ratio=1;     %  1, 1.5,   1

% range of frequencies for plant stability boundaries
Ommin=0;
Ommax=2*pi;
dOm=0.1;
Om=Ommin:dOm:Ommax;
% range of frequencies for string stability boundaries
ommin=0;
ommax=2*pi;
dom=0.00001;
om=ommin:dom:ommax;

%% Plant stability boundaries
% Omega=0
betaP1=[betamin,betamax];
alphaP1=[0,0];
% Omega>0
alphaP2=Om.^2.*cos(Om*tau)/kappa;
betaP2=Om.*sin(Om*tau)-alphaP2;

%% Sring stability boundaries
% omega=0
betaS1=[betamin,betamax];
alphaS1=[0,0];
betaS2=[betamin,betamax];
alphaS2=2*(kappa-betaS2);
% omega>0
a=(om*(kappa*tau-1)+kappa*sin(om*tau).*cos(om*tau))./...
  ((2*kappa*tau-1)*sin(om*tau)-om*tau.*cos(om*tau));
b=(om.^2.*(sin(om*tau)-om*tau.*cos(om*tau)))./...
  ((2*kappa*tau-1)*sin(om*tau)-om*tau.*cos(om*tau));
alphaS3=a+sqrt(a.^2+b);
alphaS3(a.^2+b<0)=nan;
betaS3=(om+alphaS3*kappa*tau.*sin(om*tau))./...
       (sin(om*tau)+om*tau.*cos(om*tau))-alphaS3;
betaS3(betaS3>betamax | betaS3<betamin)=nan; % remove asymptote
alphaS4=a-sqrt(a.^2+b);
alphaS4(a.^2+b<0)=nan;
betaS4=(om+alphaS4*kappa*tau.*sin(om*tau))./...
       (sin(om*tau)+om*tau.*cos(om*tau))-alphaS4;
betaS4(betaS4>betamax | betaS4<betamin)=nan; % remove asymptote

%% Stability chart
figure(1); clf; hold on; box on;
% plant stability boundaries
plot(betaP1,alphaP1,'r',betaP2,alphaP2,'r');
% string stability boundaries
plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
plot(betaS3,alphaS3,'b',betaS4,alphaS4,'b');
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with time delay',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \tau = ',num2str(tau,'%3.2f'),' [s]']);

%% Plot points on existing stability chart
% hold on; plot([-0.5,0.5,0.5],[0.4,0.4,-0.1],'bx');
% hold on; plot([0.25,0.4,0.5,0.75143,0.8],[0.4,0.4,0.4,0.4,0.4],'bx');
% axis([-0.6,1,-0.4,1.2]);