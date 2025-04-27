%% Create analytical stability charts for the vehicle chain with time delay
clear; %close all; clc;

%% Parameters
kappa=0.6;
tau0=0.6;           % delay for plant stability chart
tau=[0.6,0.65,0.7]; % delays for string stability charts
delta=-0.8:0.1:0;   % real parts of rightmost eigenvalues for plant stability contours

% range of parameters
betamin=0;         % -2,  -2,   0
betamax=1.2;          %  3,   4, 1.2
alphamin=0;        % -1,   0,   0
alphamax=1.2;         %  4,   4, 1.2
aspect_ratio=1;     %  1, 1.5,   1

% range of frequencies for plant stability boundaries
Ommin=0;
Ommax=2*pi;
dOm=0.01;
Om=(Ommin:dOm:Ommax).';
% range of frequencies for string stability boundaries
ommin=0;
ommax=2*pi;
dom=0.00001;
om=(ommin:dom:ommax).';

%% Plant stability contours with Re(s)=-delta
% Omega=0
betaP1=[betamin;betamax]+0*delta;
alphaP1=-delta./(delta+kappa).*(betaP1+delta.*exp(delta*tau0));
idx=abs(delta+kappa)<1e-6;	% delta=-kappa special case
if any(idx)
    betaP1(:,idx)=-delta(idx).*exp(delta(idx)*tau0)*[1;1];
    alphaP1(:,idx)=[alphamin;alphamax];
end
% Omega>0
alphaP2=(delta.^2+Om.^2)./(kappa*Om).*exp(delta*tau0).*...
            (delta.*sin(Om*tau0)+Om.*cos(Om*tau0));
betaP2=-exp(delta*tau0).*((delta.^2-Om.^2)./Om.*sin(Om*tau0)+...
            2*delta.*cos(Om*tau0))-alphaP2;

%% Sring stability boundaries
% omega=0
betaS1=[betamin;betamax];
alphaS1=[0;0];
betaS2=[betamin;betamax];
alphaS2=2*(kappa-betaS2);
% omega>0
a=(om.*(kappa*tau-1)+kappa*sin(om.*tau).*cos(om.*tau))./...
  ((2*kappa*tau-1).*sin(om.*tau)-om.*tau.*cos(om.*tau));
b=(om.^2.*(sin(om.*tau)-om.*tau.*cos(om.*tau)))./...
  ((2*kappa*tau-1).*sin(om.*tau)-om.*tau.*cos(om.*tau));
alphaS3=a+sqrt(a.^2+b);
alphaS3(a.^2+b<0)=nan;
betaS3=(om+alphaS3*kappa.*tau.*sin(om.*tau))./...
       (sin(om.*tau)+om.*tau.*cos(om.*tau))-alphaS3;
betaS3(betaS3>betamax | betaS3<betamin)=nan; % remove asymptote
alphaS4=a-sqrt(a.^2+b);
alphaS4(a.^2+b<0)=nan;
betaS4=(om+alphaS4*kappa.*tau.*sin(om.*tau))./...
       (sin(om.*tau)+om.*tau.*cos(om.*tau))-alphaS4;
betaS4(betaS4>betamax | betaS4<betamin)=nan; % remove asymptote

%% Stability chart
% plant stability boundaries
figure(11); clf; hold on; box on;
plot(betaP1,alphaP1,'r',betaP2,alphaP2,'r');
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with time delay',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \tau = ',num2str(tau0,'%3.2f'),...
       ' [s],   \delta = [',regexprep(num2str(delta),'\s+',', '),']']);

% string stability boundaries
figure(12); clf; hold on; box on;
plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
plot(betaS3,alphaS3,'b',betaS4,alphaS4,'b');
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with time delay',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \tau = [',regexprep(num2str(tau),'\s+',', '),'] [s]']);