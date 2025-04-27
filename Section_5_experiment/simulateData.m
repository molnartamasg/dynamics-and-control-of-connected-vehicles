%% Perform simulations for a 2 vehicle chain using data as input
clear; %close all; clc;
addpath('data');

%% Select data set to be loaded
% 2 vehicles, 1 CHV + 1 CAV
% #1 - a=0.4, b1=0.5
load('2vehicles_1.mat');	t0=6;	kCAV=2;	dh=0;	alpha=0.4;	beta=0.5;
tminplot=0; tmaxplot=150; vminplot=0; vmaxplot=15; aminplot=-4; amaxplot=4; hminplot=0; hmaxplot=35;
% #2 - a=0.6, b1=0.8
% load('2vehicles_2.mat');	t0=8.9;	kCAV=2; dh=0;   alpha=0.6;	beta=0.8;
% tminplot=0; tmaxplot=150; vminplot=0; vmaxplot=15; aminplot=-4; amaxplot=4; hminplot=0; hmaxplot=35;

% 8 vehicles, 6 CHV + 1 CAV + 1 CHV
% #1 - a=0.4, b1=0.5, b2=0, b3=0
% load('8vehicles_1.mat');	t0=270;	kCAV=7;	dh=3;	alpha=0.4;	beta=[0.5,0,0,0,0,0];
% tminplot=0; tmaxplot=500; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;
% #4 - a=0.4, b1=0.2, b2=0.3, b3=0.3
% load('8vehicles_4.mat');	t0=330;	kCAV=7;	dh=3;	alpha=0.4;	beta=[0.2,0.3,0.3,0,0,0];
% tminplot=0; tmaxplot=500; vminplot=0; vmaxplot=30; aminplot=-10; amaxplot=5; hminplot=0; hmaxplot=100;

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

%% Parameters for simulation
% fixed parameters for CAV
sigma=0.6;  % 0.6, 0.7, 0.8
hst=5;
hgo=55;
vmax=30;
amin=7;
amax=3;
betasum=sum(beta);

% correct headway to the one the vehicle perceived
% dh=0;	% for safety reasons, the actual controller of the CAV
        % assumed a larger vehicle length than the actual one (indicated by dh)
        % and calculated a smaller headway accordingly
        % so the perceived headway was smaller than the real one in the data
        % and we will use the perceived one
hdwy{kCAV}=hdwy{kCAV}-dh;
        
% constants for loss terms and vehicle dynamics
% gamma=0.01;         % [-] tyre rolling resistance coefficient
% g=9.81;             % [m/s^2] gravitatioinal constant
% a=gamma*g;          % [m/s^2]
% k=0.45;             % air drag [kg/m]
m=1500;             % mass [kg]
% c=k/m;              % [1/m]
P=75000;            % [W] max power
% P=50000;            % [W] max power

% range policy, speed policy, and saturations for CAV
V=@(h)vmax*(hgo<=h) + vmax*(h-hst)/(hgo-hst).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
% sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
amaxv=@(v)min([amax+0*v,P/m./abs(v)],[],2);
sat=@(v,u)(u<-amin).*(-amin)+(-amin<=u & u<=amaxv(v)).*u+(amaxv(v)<u).*amaxv(v);

% control input for CAV
u=@(h,v,vL)alpha*V(h)-(alpha+betasum)*v+beta*W(vL);

% initial conditions for CAV
h0=interp1(tCom{kCAV},hdwy{kCAV},t0);
v0=interp1(time{kCAV},vel{kCAV},t0);
xinit=@(t)[h0;v0];

% simulation time
tend=max(vertcat(time{:}));
deltat=min(diff(time{1}));
tsim=(t0:deltat:tend).';        

% number of preceding CHVs that the CAV may respond to
chv_num=length(beta);
% preceding CHVs' velocity from data
vLead=nan(length(tsim),chv_num);
vLeaddelay=nan(size(vLead));
for kL=1:chv_num
    vLead(:,kL)=interp1(time{kCAV-kL},vel{kCAV-kL},tsim,'linear','extrap');
    vLeaddelay(:,kL)=interp1(tsim,vLead(:,kL),tsim-sigma,'linear','extrap');
end
vL=@(t)vLead(t==tsim,:).';
vLdelay=@(t)vLeaddelay(t==tsim,:).';
select=@(vL,kL)vL(kL);  % function for selecting element from vector

% title to put on figure
problem='CCC using measurement data';
% list of parameters to put on figure
parlist=['\alpha=',num2str(alpha,'%3.2f'),' [1/s]   ',...
         '\beta=[',regexprep(num2str(beta),'\s+',', '),'] [1/s]   '...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];
% legend to put on figure
vehtype=repmat({'measured CHV'},1,veh_num);
vehtype{kCAV}='measured CAV';
vehlegend=horzcat(vehtype,{'simulated CAV'});

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay)[select(vL(t),1)-x(2);
                    -a-c*x(2)^2+sat(x(2),u(xdelay(1),xdelay(2),vLdelay(t)))];

% perform simulation
x=ddeab4(@(t,x,xdelay)model(t,x,xdelay),sigma,xinit,tsim);

% extract headway and velocity
hdwysim=x(1,:).';
velsim=x(2,:).';

% calculate accelearation
accsim=sgolayfilt(diff(velsim)/deltat,3,21);
accsim=[accsim(1);accsim];
    
% calculate energy consumption
wsim=cumsum(max(accsim+a+c*velsim.^2,0).*velsim*deltat);

%% Plot of solution
% plot velocity vs time
figure(1); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured velocity
for kk=1:veh_num
   LL(kk)=plot(time{kk}-t0,vel{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
% plot speed of simulated CAV
LL(end)=plot(tsim-t0,velsim,'Linewidth',2.5,'Color','r');
xlim([tminplot,tmaxplot]);
ylim([vminplot,vmaxplot]);
xlabel('Time [s]');ylabel('Speed [m/s]');
title(['Simulation results of ',problem,10,parlist]);
legend(LL,vehlegend,'Location','northwest');

% plot acceleration vs time
figure(2); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured acceleration
for kk=1:veh_num
   LL(kk)=plot(time{kk}-t0,acc{kk},'Linewidth',1.5,'Color',colours(kk,:));
%    LL(kk)=plot(ttt{kk}-t0,aaa{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
% plot acceleration of simulated vehicle
LL(end)=plot(tsim-t0,accsim,'Linewidth',2.5,'Color','r');
xlim([tminplot,tmaxplot]);
ylim([aminplot,amaxplot]);
xlabel('Time [s]');ylabel('Acceleration [m/s^2]');
title(['Simulation results of ',problem,10,parlist]);
legend(LL,vehlegend,'Location','southwest');

% plot headway vs time
figure(3); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured headway
for kk=1:veh_num
   LL(kk)=plot(tCom{kk}-t0,hdwy{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
% plot headway of simulated vehicle
LL(end)=plot(tsim-t0,hdwysim,'Linewidth',2.5,'Color','r');
xlim([tminplot,tmaxplot]);
ylim([hminplot,hmaxplot]);
xlabel('Time [s]');ylabel('Headway [m]');
title(['Simulation results of ',problem,10,parlist]);
legend(LL,vehlegend,'Location','northwest');

% plot energy consumption vs time
figure(4); clf; hold on; box on;
LL=zeros(veh_num+1,1);
% plot measured energy consumption
for kk=1:veh_num
   LL(kk)=plot(time{kk}-t0,encons{kk},'Linewidth',1.5,'Color',colours(kk,:));
%    LL(kk)=plot(ttt{kk}-t0,www{kk},'Linewidth',1.5,'Color',colours(kk,:));
end
% plot energy consumption of simulated vehicle
LL(end)=plot(tsim-t0,wsim,'Linewidth',2.5,'Color','r');
xlim([tminplot,tmaxplot]);
xlabel('Time [s]');ylabel('Energy consumption [J/kg]');
title(['Simulation results of ',problem,10,parlist]);
legend(LL,vehlegend,'Location','northwest');

%% Evaluation of frequency response
% for FFT, consider only those times when vehicles were moving
% based on the specified start time t0
% get time interval (t1<=t<=t2) corresponding to moving
t1=t0;
t2=time{kCAV}(end);

% pick time vector to which data is interpolated before doing FFT
tFFT=tsim(t1<=tsim & tsim<=t2);
% get number of data points
n=length(tFFT);
% get sampling frequency and discrete frequency values
fs=1/deltat;
freq=(0:floor(n/2))/n*fs;

% frequency content of measured velocity
velFFT=cell(size(vel));
for kk=1:veh_num
    % interpolate velocities to the same time vector
    velinterp=interp1(time{kk},vel{kk},tFFT,'linear','extrap');
    % calculate FFT for velocity fluctuations
    velFFT{kk}=abs(fft(velinterp-mean(velinterp))/n);
    % take twice the first half of the FFT
    velFFT{kk}=velFFT{kk}(1:floor(n/2+1));
    velFFT{kk}(2:end)=2*velFFT{kk}(2:end);
    % filter the FFT for smoother results
    velFFT{kk}=abs(sgolayfilt(velFFT{kk},3,31));
end

% frequency content of simulated velocity
velinterp=interp1(tsim,velsim,tFFT,'linear','extrap');
velsimFFT=abs(fft(velinterp-mean(velinterp))/n);
velsimFFT=velsimFFT(1:floor(n/2+1));
velsimFFT(2:end)=2*velsimFFT(2:end);
velsimFFT=abs(sgolayfilt(velsimFFT,3,31));

% frequency range for plotting and string instability index
fmin=0;
fmax=1;

% plot frequency content
figure(5); clf; box on;
semilogy([fmin,fmax],[1,1],'k--','HandleVisibility','off'); hold on;
for kk=1:veh_num
    semilogy(freq,velFFT{kk},'Linewidth',1.5,'Color',colours(kk,:));
    hold on;
end
semilogy(freq,velsimFFT,'Linewidth',2.5,'Color','r');
xlim([fmin,fmax]);
ylim([0.0001,10]);
xlabel('Frequency [Hz]');ylabel('Velocity [m/s]');
title(['Frequency content of velocity fluctuations',10,parlist]);
legend(vehlegend,'Location','northeast');

% plot frequency response wrt the head vehicle of the full group
figure(6); clf; box on;
semilogy([fmin,fmax],[1,1],'k--','HandleVisibility','off'); hold on;
for kk=2:veh_num
    semilogy(freq,velFFT{kk}./velFFT{1},'Linewidth',1.5,'Color',colours(kk,:));
    hold on;
end
semilogy(freq,velsimFFT./velFFT{1},'Linewidth',2.5,'Color','r');
xlim([fmin,fmax]);
ylim([0.1,10]);
xlabel('Frequency [Hz]');ylabel('Amplification');
title(['Frequency repsonse w.r.t. the head vehicle of the full group',10,parlist]);
legend(vehlegend(2:end),'Location','northeast');

% string stability index
fmini = find(freq>=fmin,1);
fmaxi = find(freq>=fmax,1);
freqidx=fmini:fmaxi;
Cs=sum(max(velFFT{kCAV}(freqidx)./velFFT{1}(freqidx)-1,0)*fs/n)/(fmax-fmin);
% Cs=sum(max(velsimFFT(freqidx)./velFFT{1}(freqidx)-1,0)*fs/n)/(fmax-fmin);