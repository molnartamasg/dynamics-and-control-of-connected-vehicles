%% Plot experimental data
clear; %close all; clc;
addpath('data');

% Main data:
% time - time
% arc - position (arclength, starts from zero for the leading car)
% vel - velocity
% acc - acceleration
% encons - energy consumption
% tCom - common time instants when data is available for both vehicles in a pair
% hdwy - headway from predecessor vehicle at these common time instants
%        (great circle distance minus 5 meters vehicle length)

%% Select data set to be loaded
% 2 vehicles, 1 CHV + 1 CAV
% #1 - a=0.4, b1=0.5
% load('2vehicles_1.mat');	t0=6;	kCAV=2;
% tminplot=0; tmaxplot=150; vminplot=0; vmaxplot=15; aminplot=-4; amaxplot=4; hminplot=0; hmaxplot=35;
% #2 - a=0.6, b1=0.8
% load('2vehicles_2.mat');	t0=8.9;	kCAV=2;
% tminplot=0; tmaxplot=150; vminplot=0; vmaxplot=15; aminplot=-4; amaxplot=4; hminplot=0; hmaxplot=35;

% 2 vehicles, 1 CHV + 1 CAV
% #3 - human
% load('2vehicles_3.mat');	t0=1;
% tminplot=0; tmaxplot=25; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=160;
% #4 - a=0.4, b1=0.5
% load('2vehicles_4.mat');	t0=4;	kCAV=2;	dh=3;
% tminplot=0; tmaxplot=25; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=160;

% 3 vehicles, 2 CHV + 1 CAV
% #1 - human
% load('3vehicles_1.mat');	t0=40;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #2 - a=0.4, b1=0.5, b2=0
% load('3vehicles_2.mat');	t0=43;	kCAV=3;	dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #3 - a=0.4, b1=0.2, b2=0.6
% load('3vehicles_3.mat');	t0=79;	kCAV=3;	dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;

% 4 vehicles, 3 CHV + 1 CAV
% #1 - human
% load('4vehicles_1.mat');	t0=32;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #2 - a=0.4, b1=0.5, b2=0, b3=0
% load('4vehicles_2.mat');	t0=42;  kCAV=4; dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #3 - a=0.4, b1=0.2, b2=0.3, b3=0
% load('4vehicles_3.mat');	t0=46;  kCAV=4; dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #4 - a=0.4, b1=0.2, b2=0, b3=0.3
% load('4vehicles_4.mat');	t0=41;  kCAV=4; dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;
% #5 - a=0.4, b1=0.2, b2=0.3, b3=0.3
% load('4vehicles_5.mat');	t0=44;  kCAV=4; dh=3;
% tminplot=0; tmaxplot=10; vminplot=0; vmaxplot=20; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=50;

% 8 vehicles, 6 CHV + 1 CAV + 1 CHV
% #1 - human
% load('8vehicles_1.mat');	t0=0;
% tminplot=60; tmaxplot=560; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;	% 60, 315, 355, 560
% #2 - a=0.4, b1=0.5, b2=0, b3=0
% load('8vehicles_2.mat');	t0=0;	kCAV=7;	dh=3;
% tminplot=270; tmaxplot=770; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;	% 270, 375, 415, 770
% #3 - a=0.4, b1=0.2, b2=0.3, b3=0
% load('8vehicles_3.mat');	t0=0;	kCAV=7;	dh=3;
% tminplot=300; tmaxplot=800; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;	% 300, 450, 490, 800
% #4 - a=0.4, b1=0.2, b2=0, b3=0.3
% load('8vehicles_4.mat');	t0=0;	kCAV=7;	dh=3;
% tminplot=30; tmaxplot=530; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;	% 30, 375, 415, 530
% #5 - a=0.4, b1=0.2, b2=0.3, b3=0.3
% load('8vehicles_5.mat');	t0=0;	kCAV=7;	dh=3;
% tminplot=330; tmaxplot=830; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;	% 330, 335, 375, 830

% 3 vehicles, 1 CAV + 2 CHV (virtual ring)
% #1 - a=0.4, b1=0.5, b2=0, k=0.6
% load('virtualring_1.mat');	t0=45;	kCAV=1;
% tminplot=0; tmaxplot=100; vminplot=5; vmaxplot=25; aminplot=-6; amaxplot=4; hminplot=0; hmaxplot=60;
% #2 - a=0.4, b1=0.5, b2=0, k=0.8
% load('virtualring_2.mat');	t0=90;	kCAV=1;
% tminplot=0; tmaxplot=100; vminplot=5; vmaxplot=25; aminplot=-6; amaxplot=4; hminplot=0; hmaxplot=60;
% #3 - a=0.4, b1=0.5, b2=0, k=1.0
% load('virtualring_3.mat');	t0=85;	kCAV=1;
% tminplot=0; tmaxplot=100; vminplot=5; vmaxplot=25; aminplot=-6; amaxplot=4; hminplot=0; hmaxplot=60;
% #4 - a=0.4, b1=0.2, b2=0.3, k=1.0
load('virtualring_4.mat');	t0=20;	kCAV=1;
tminplot=0; tmaxplot=100; vminplot=5; vmaxplot=25; aminplot=-6; amaxplot=4; hminplot=0; hmaxplot=60;

%% Process data
veh_num=size(colours,1);

% Compute acceleration with interpolation
dt=0.1;
ttt=cell(size(time));
vvv=cell(size(vel));
aaa=cell(size(acc));
for kk=1:veh_num
    ttt{kk}=(round(time{kk}(1),1):dt:round(time{kk}(end),1)).';
    vvv{kk}=interp1(time{kk},vel{kk},ttt{kk},'linear','extrap');
    aaa0=sgolayfilt(diff(vvv{kk})/dt,3,21);
    aaa{kk}=[aaa0(1);aaa0];
end

% Compute energy consumption according to acceleration
a=9.81e-2;
c=3e-4; % 3e-4, 2.74e-4
www=cell(size(encons));
for kk=1:veh_num
    www{kk} = cumsum(max(aaa{kk}+a+c*vvv{kk}.^2,0).*vvv{kk}*dt);
end

%% Plot data
% Plot position vs time
figure(1);clf;hold on;box on;grid on;
for kk=1:veh_num
    plot(time{kk}-t0,arc{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
xlabel('Time [s]');ylabel('Position [m]');
xlim([tminplot,tmaxplot]);

% Plot velocity vs time
figure(2);clf;hold on;box on;grid on;
for kk=1:veh_num
    plot(time{kk}-t0,vel{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
xlabel('Time [s]');ylabel('Speed [m/s]');
xlim([tminplot,tmaxplot]);
ylim([vminplot,vmaxplot]);

% Plot acceleration vs time
figure(3);clf;hold on;box on;grid on;
for kk=1:veh_num
    plot([tminplot,tmaxplot],[0,0],'k--','Linewidth',1);
%     plot(time{kk}-t0,acc{kk},'Linewidth',1.5,'Color',colours(kk,:));
    plot(ttt{kk}-t0,aaa{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
xlim([tminplot,tmaxplot]);
ylim([aminplot,amaxplot]);

% Plot headway vs time
figure(4);clf;hold on;box on;grid on;
for kk=1:veh_num
    % adjust headway according to vehicle length buffer
    if exist('kCAV','var') && exist('dh','var')&& kk==kCAV
        plot(tCom{kk}-t0,hdwy{kk}-dh,'Linewidth',1.5,'Color',colours(kk,:));
    else
        plot(tCom{kk}-t0,hdwy{kk},'Linewidth',1.5,'Color',colours(kk,:));
    end
end
xlabel('Time [s]');ylabel('Headway [m]');
xlim([tminplot,tmaxplot]);
ylim([hminplot,hmaxplot]);

% Plot energy consumption vs time
figure(5);clf;hold on;box on;grid on;
for kk=1:veh_num
%     plot(time{kk}-t0,encons{kk},'Linewidth',1.5,'Color',colours(kk,:));
    plot(ttt{kk}-t0,www{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
xlabel('Time [s]');ylabel('Energy consumption [J/kg]');
xlim([tminplot,tmaxplot]);