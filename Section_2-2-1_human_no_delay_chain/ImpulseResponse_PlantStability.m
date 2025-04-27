%% Calculate impulse response of the car following model

clear; clc; % close all;

%% Set parameters
% parameters of car following model
alpha=2;  % 2 -0.1
beta=0.0;
kappa=0.6;

%% Impulse response
% system matrices
a=[0,-1;alpha*kappa,-(alpha+beta)];
b=[1;beta];
c=[0,1];
A=1;

% time
t0=0;
tend=5;
deltat=0.01;
time=t0:deltat:tend;

% symbolic response
syms t tau
cexpatb=c*expm(a*t)*b;
y1sym=A*cexpatb;
% y2sym=int(subs(cexpatb,t,t-tau)*subs(y1sym,t,tau),tau,0,t);
% y3sym=int(subs(cexpatb,t,t-tau)*subs(y2sym,t,tau),tau,0,t);
% y4sym=int(subs(cexpatb,t,t-tau)*subs(y3sym,t,tau),tau,0,t);

% numerical response
y1num=eval(subs(y1sym,t,time));
% y2num=eval(subs(y2sym,t,time));
% y3num=eval(subs(y3sym,t,time));
% y4num=eval(subs(y4sym,t,time));

%% Plot response
figure(1); clf; hold on; box on;
plot(time,y1num);
% plot(time,y2num);
% plot(time,y3num);
% plot(time,y4num);
xlabel('Time [s]');
ylabel('Response [m/s]');
title(['Impulse response of vehicle chain assuming OVM without delay',10,...
        '\alpha = ',num2str(alpha,'%3.2f'),...
        ' [1/s],   \beta = ',num2str(beta,'%3.2f'),' [1/s]']);
axis([t0,tend,-A*0.5,A*0.5]);