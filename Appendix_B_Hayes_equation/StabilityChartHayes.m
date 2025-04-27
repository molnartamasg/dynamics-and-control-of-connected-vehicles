%% Create analytical stability charts for the Hayes equation
clear; %close all; clc;

%% Parameters
tau=1;
% range of parameters
amin=-15;
amax=15;
bmin=-15;
bmax=15;
aspect_ratio=1;

% range of frequencies for stability boundaries
ommin=0;
ommax=10*pi;
dom=0.001;
om=ommin:dom:ommax;

%% Stability boundaries
% omega=0
a1=[amin,amax];
b1=-a1;
% omega>0
a2=om./tan(om*tau);
b2=-om./sin(om*tau);
% remove asymptotes
b2(b2>bmax | b2<bmin)=nan;

%% Stability chart
figure(1); clf; hold on; box on;
plot([amin,amax],[0,0],'k--');
plot([0,0],[bmin,bmax],'k--');
plot(a1,b1,'r',a2,b2,'r');
axis([amin amax bmin bmax]);
pbaspect([aspect_ratio,1,1]);
xlabel('a');
ylabel('b');
title(['Stability chart of the Hayes equation with',...
       '\tau = ',num2str(tau,'%3.2f')]);