%% Create analytical stability charts for ACC
clear; %close all; clc;

%% Parameters
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
% d=0;

% control parameters (range policy and integral gain)
eta=0.1;	% 0.0, 0.05, 0.1, 0.2, 0.8
kappa=0.6;  % 0.4, 0.6, 0.8, 1.0
sigma=0.6;  % 0.0, 0.6

% range of parameters
betamin=-2;     % -2,   0
betamax=3;      %  3, 1.2
alphamin=-1;	% -1,   0
alphamax=4;     %  4, 1.2
aspect_ratio=1;	%  1,   1

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
% % expression of the s=0 plant stability boundary
% if there is I term
% eta=0;
% if there is no I term
if eta==0
    betaP1=[betamin,betamax];
    alphaP1=[0,0];
end

% expression of the s=jOm plant stability boundary
alphaP2=(Om.^2.*cos(Om*sigma)+d*Om.*sin(Om*sigma)-eta)/kappa;
betaP2=Om.*sin(Om*sigma)-d*cos(Om*sigma)+eta*kappa./Om.^2-alphaP2;

% % expressions of the om=0 string stability boundaries
% if there is I term
% eta=0;
% eta=2*d*kappa;
% if there is no I term
if eta==0
    betaS1=linspace(betamin,betamax,10001);
    alphaS1=-(betaS1+d-kappa-kappa*d*sigma)-...
                sqrt((betaS1+d-kappa-kappa*d*sigma).^2-d*(d+2*betaS1));
    alphaS1(imag(alphaS1)~=0)=nan;
    betaS2=linspace(betamin,betamax,10001);
    alphaS2=-(betaS2+d-kappa-kappa*d*sigma)+...
                sqrt((betaS2+d-kappa-kappa*d*sigma).^2-d*(d+2*betaS2));
    alphaS2(imag(alphaS2)~=0)=nan;
end

% expressions of om>0 string stability boundaries
if sigma>0
    c0=om.^4-eta^2+eta*((kappa-d)*om.^2*sigma+2*kappa*d).*cos(om*sigma)...
        +eta*(om.^2*sigma-(kappa-d)+kappa*d*sigma).*om.*sin(om*sigma);
    c1=kappa*(om.^2*sigma+d).*om.*sin(om*sigma)-kappa*d*sigma*om.^2.*cos(om*sigma);
    c2=(d*sigma+1)*om.^3.*sin(om*sigma)+om.^4*sigma.*cos(om*sigma);
    aa=(-c0+(c2*kappa-c1*d).*cos(om*sigma)...
        +(c1+kappa*d*c2./om.^2).*om.*sin(om*sigma))./(2*c1-c2);
    bb=(-c2.*(om.^2+d^2+eta^2./om.^2)+2*(c2*eta-c0*d+c2*eta*kappa*d./om.^2).*cos(om*sigma)...
        +2*(c0-c2*eta*(kappa-d)./om.^2).*om.*sin(om*sigma))./(2*c1-c2);
    alphaS3=aa+sqrt(aa.^2+bb);
    alphaS3(aa.^2+bb<0)=nan;
    betaS3=(c0+c1.*alphaS3)./c2-alphaS3;
    betaS3(betaS3>betamax | betaS3<betamin)=nan; % remove asymptote
    alphaS4=aa-sqrt(aa.^2+bb);
    alphaS4(aa.^2+bb<0)=nan;
    betaS4=(c0+c1.*alphaS4)./c2-alphaS4;
    betaS4(betaS4>betamax | betaS4<betamin)=nan; % remove asymptote
else % treat zero delay case separately
    betaS3=linspace(betamin,betamax,10001);
    alphaS3=-(d+betaS3-kappa)+sqrt((betaS3-kappa).^2-2*d*kappa+2*eta-2*sqrt(eta*(eta-2*d*kappa)));
    alphaS3(imag(alphaS3)~=0)=nan;
    betaS4=linspace(betamin,betamax,10001);
    alphaS4=-(d+betaS4-kappa)-sqrt((betaS4-kappa).^2-2*d*kappa+2*eta-2*sqrt(eta*(eta-2*d*kappa)));
    alphaS4(imag(alphaS4)~=0)=nan;
end

%% Stability chart
figure(1); clf; hold on; box on;
% plot s=0 plant stability boundary
if eta==0
    plot(betaP1,alphaP1,'r');
end
% plot s=jOm plant stability boundary
plot(betaP2,alphaP2,'r');
% plot om=0 string stability boundaries
if eta==0
    plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
end
% plot om>0 string stability boundaries
plot(betaS3,alphaS3,'b',betaS4,alphaS4,'b');
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\alpha [1/s]');
ylabel('\beta [1/s]');
title(['Stability chart of ACC',10,...
       vehtype,'   \eta=',num2str(eta,'%3.2f'),...
       ' [1/s^2],   \kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \sigma = ',num2str(sigma,'%3.2f'),' [s]']);