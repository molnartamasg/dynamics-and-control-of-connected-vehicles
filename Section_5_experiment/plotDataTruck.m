%% Plot experimental data
clear; %close all; clc;
addpath('data');

% Main data:
% t - time
% v - truck's velocity
% v1 - lead vehicle's velocity
% a - truck's acceleration
% a1 - lead vehicle's acceleration
% h - headway

%% Select data set to be loaded
% 2 vehicles, 1 CHV + 1 CAV (truck)
% #1 - kappa=0.8, alpha=0.4, beta=0.5, gamma=0
% load('truck_1.mat');
% tminplot=184; tmaxplot=194; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=30;
% #2 - kappa=0.8, alpha=0.4, beta=0.5, gamma=0.25
load('truck_2.mat');
tminplot=184; tmaxplot=194; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=30;

%% Plot data
% Plot velocity vs time
figure(2);clf;hold on;box on;grid on;
plot(t,v1,'Linewidth',1.5,'Color',[1,0,0]);
plot(t,v,'Linewidth',1.5,'Color',[0.5,0,0.5]);
xlabel('Time [s]');ylabel('Speed [m/s]');
xlim([tminplot,tmaxplot]);
ylim([vminplot,vmaxplot]);

% Plot acceleration vs time
figure(3);clf;hold on;box on;grid on;
plot(t,a1,'Linewidth',1.5,'Color',[1,0,0]);
plot(t,a,'Linewidth',1.5,'Color',[0.5,0,0.5]);
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
xlim([tminplot,tmaxplot]);
ylim([aminplot,amaxplot]);

% Plot headway vs time
figure(4);clf;hold on;box on;grid on;
plot(t,h,'Linewidth',1.5,'Color',[0.5,0,0.5]);
xlabel('Time [s]');ylabel('Headway [m]');
xlim([tminplot,tmaxplot]);
ylim([hminplot,hmaxplot]);