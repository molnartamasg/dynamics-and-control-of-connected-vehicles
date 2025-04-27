%% Create analytical stability charts for the vehicle chain with ZOH
clear; %close all; clc;

%% Parameters
kappa=0.6;
deltat=0.4;	% 0.1, 0.4

% range of parameters
betamin=-2;
betamax=3;
alphamin=-1;
alphamax=4;
aspect_ratio=1;

% range of frequencies for plant stability boundaries
Ommin=-pi/deltat;
Ommax=pi/deltat;
dOm=2*pi/deltat/200;
Om=Ommin:dOm:Ommax;

% range of frequencies for string stability boundaries
ommin=0;
ommax=2*pi/deltat;
dom=2*pi/deltat/10000;

%% Plant stability boundaries
% z=1
betaP1=[betamin,betamax];
alphaP1=[0,0];
% z=-1
betaP2=[betamin,betamax];
alphaP2=-2/deltat-betaP2;
% z=exp(j*Om*dt)
alphaP3=2/kappa/deltat^2*(2*sin(Om*deltat).^2+3*cos(Om*deltat)-3);
betaP3=1/deltat*(2*sin(Om*deltat).^2+cos(Om*deltat)-1)-alphaP3;

%% Sring stability boundaries
% om=0
betaS1=[betamin,betamax];
alphaS1=[0,0];
betaS2=[betamin,betamax];
alphaS2=2*(kappa-betaS2);
% om=pi/dt
betaS3=[betamin,betamax];
alphaS3=-2/deltat+0*betaS3;
betaS4=[betamin,betamax];
alphaS4=-2/deltat-2*betaS4;
% 0<om~=pi/dt
syms k dt a b z om real
T=((dt^2/2*a*k+dt*b)*z+dt^2/2*a*k-dt*b)/...
    (z^3-2*z^2+(dt^2/2*a*k+dt*(a+b)+1)*z+dt^2/2*a*k-dt*(a+b));
Tom=subs(T,z,exp(1i*om*dt));
[N,D]=numden(T);
[Nom,Dom]=numden(Tom);
P=simplify(real(Dom)^2+imag(Dom)^2-real(Nom)^2-imag(Nom)^2)/om^2;
Q=diff(P,om)/2/om;
% get symbolic solutions for the stability boundaries in the form
% alpha=alpha(om), beta=beta(om)
solu=solve([P==0,Q==0],[b,a],'MaxDegree',3);

om=ommin:dom:ommax;

%% Stability chart
figure(1); clf; hold on; box on;
% plant stability boundaries (z=1, z=-1, z=exp(j*Om*dt))
plot(betaP1,alphaP1,'r',betaP2,alphaP2,'r',betaP3,alphaP3,'r');
% om=0 string stability boundaries
plot(betaS1,alphaS1,'b',betaS2,alphaS2,'b');
% om=pi/dt string stability boundaries
plot(betaS3,alphaS3,'b',betaS4,alphaS4,'b');
% om>0 string stability boundaries
for ksolu=1:length(solu.b)
    % evaluate the string stability limit for different values of omega
    betastring=eval(subs(solu.b(ksolu),[k,dt],[kappa,deltat]));
    alphastring=eval(subs(solu.a(ksolu),[k,dt],[kappa,deltat]));
    % remove complex values
    betastring(imag(betastring)~=0)=nan;
    alphastring(imag(alphastring)~=0)=nan;
    % remove possible asymptotes
    betastring(betastring>betamax | betastring<betamin)=nan;
    % plot omega>0 string stability boundaries
    plot(betastring,alphastring,'b');
end
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the vehicle chain with ZOH',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \Deltat = ',num2str(deltat,'%3.2f'),' [s]']);

%% Plot points on existing stability chart
% deltat=0.4 case
% hold on; plot([-0.5,0.25,0.5,0.8,0.5],[0.4,0.4,0.4,0.4,-0.1],'bx');