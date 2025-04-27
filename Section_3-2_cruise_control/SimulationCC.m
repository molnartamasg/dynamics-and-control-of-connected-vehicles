%% Perform simulations for cruise controller
clear; %close all; clc;

%% Parameters
% type of vehicle
vehtype='car';
% vehtype='truck';

% vehicle parameters
if strcmp(vehtype,'car')	% personal vehicle
    k=0.45;     % air drag
    m=1500;     % mass
    amin=7;     % max deceleration
    amax=3;     % max acceleration
    P=75000;	% max power
else	% truck
    k=3.85;     % air drag
    m=29500;	% mass
    amin=6;     % max deceleration
    amax=2;     % max acceleration
    P=300000;	% max power
end

% gravitational constant, rolling resistance, elevation
g=9.81;
gamma=0.01;
phi=0.025;

% resistance terms
p=@(v)g*sin(phi)+gamma*g*cos(phi).*sign(v)+k/m*v.^2.*sign(v);
% acceleration constraints
a=@(v)min([amax+0*v;P/m./abs(v)]);
% saturation function
sat=@(v,u)(u<-amin).*(-amin)+(-amin<=u & u<=a(v)).*u+(a(v)<u).*a(v);

% control gains
alpha=0.9;
eta=0.8;	% 0.0, 0.05, 0.2, 0.8

% initial and reference velocity
v0=15;
vr=20;

% regions for plotting
vminplot=15;
vmaxplot=25;

% time
t0=0;
tend=20;	% 90, 20
deltat=0.01;

% range of spectrum plot
Remin=-1;
Remax=0.2;
Immin=-1;
Immax=1;

%% Calculate steady state error
vstar=(eta~=0)*vr+(eta==0)*fsolve(@(vstar)alpha*(vr-vstar)-p(vstar),vr);
ess=vr-vstar;
disp(['Steady state error: ',num2str(ess), ' [m/s]']);

%% Simulation

% equation of system and cruise control
cruise=@(t,x)[-p(x(1))+sat(x(1),alpha*(vr-x(1))+eta*x(2));...
              vr-x(1)];

% simulation time
time=t0:deltat:tend;	% simulation time
% initial state
x0=[v0;0];
% perform simulation
sol=ode45(@(t,x)cruise(t,x),[t0,tend],x0);
x=deval(sol,time);

%% Plot of solution
% calculate parameters of linearized system
zeta=(2*k/m*vr+alpha)/2/sqrt(eta);
omn=sqrt(eta);
omd=sqrt(1-zeta^2)*omn;
s=[-zeta*omn-1i*omd,-zeta*omn+1i*omd];
if eta==0
    s=[-(2*k/m*vr+alpha),0];
end

% plot related spectrum
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
plot(real(s),imag(s),'b*')
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(s)');
ylabel('Im(s)');
title(['Spectrum of the cruise controller',10,...
       vehtype,'   \alpha=',num2str(alpha),' [1/s]','   \eta=',num2str(eta),' [1/s^2]']);
   
% plot simulation results
figure(2); clf; hold on; box on;
if zeta<1   % fit exponential envelope to simulated signal
    tstar=2*pi/omd/2;
    [~,idx]=min(abs(time-tstar));
    vstar=x(1,idx);
    tstar=time(idx);
	plot(time,vr+(vstar-vr)*exp(-zeta*omn*(time-tstar)),'k--');
	plot(time,vr-(vstar-vr)*exp(-zeta*omn*(time-tstar)),'k--');
end
plot(time,x(1,:),'r');
plot([t0,tend],[vr,vr],'k--');
axis([t0 tend vminplot vmaxplot]);
xlabel('t [s]');
ylabel('v [m/s]');
title(['Simulations for the cruise controller',10,...
       vehtype,'   \alpha=',num2str(alpha),' [1/s]','   \eta=',num2str(eta),...
       ' [1/s^2]   e_{ss}=',num2str(ess,'%4.3f'),' [m/s]']);
   
% plot corresponding acceleration
figure(3); clf; hold on; box on;
plot(time,gradient(x(1,:))./gradient(time),'r');
plot([t0,tend],[0,0],'k--');
axis([t0 tend -amin amax]);
xlabel('t [s]');
ylabel('a [m/s^2]');
title(['Simulations for the cruise controller',10,...
vehtype,'   \alpha=',num2str(alpha),' [1/s]','   \eta=',num2str(eta),...
' [1/s^2]   e_{ss}=',num2str(ess,'%4.3f'),' [m/s]']);