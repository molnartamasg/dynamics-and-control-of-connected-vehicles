%% Perform simulations for CCC in a 4 vehicle scenario
clear; %close all; clc;

%% Parameters
% fixed parameters
hst=5;
hgo=55;
vmax=30;
vstar=20;
amin=7;
amax=3;
l=5;
% CHV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;
% CAV
kappa0=0.6;
alpha01=0.4;
beta01=0.5;
beta02=0.5; % 0.0, 0.5, 0.0, 0.5
beta03=0.5; % 0.0, 0.0, 0.5, 0.5
sigma=0.6;

% range policy
% CHV
% Vh=@(h)vmax*(hgo<=h)...
%     + vmax*(2*hgo-hst-h).*(h-hst)/(hgo-hst)^2.*(hst<h & h<hgo);
% hstarh=hgo-sqrt(1-vstar/vmax)*(hgo-hst);
Vh=@(h)vmax*(hgo<=h)...
    + kappah*(h-hst).*(hst<h & h<hgo);
hstarh=hst+vstar/kappah;
% CAV
V0=@(h)vmax*(hgo<=h)...
    + kappa0*(h-hst).*(hst<h & h<hgo);
hstar0=hst+vstar/kappa0;

% saturation function and speed policy
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);

% control input
% CHV
uh=@(s,v,sL,vL)alphah*(Vh(sL-s-l)-v)+betah*(vL-v);
% CAV
u0=@(s0,v0,s1,v1,v2,v3)alpha01*(V0(s1-s0-l)-v0)+...
    beta01*(W(v1)-v0)+beta02*(W(v2)-v0)+beta03*(W(v3)-v0);

% regions for plotting
sminplot=-200;
smaxplot=1200;
vminplot=8;
vmaxplot=22;
hminplot=20;
hmaxplot=45;

% title to put on figure
problem='CCC with 4 vehicles';
% list of parameters to put on figure
parlist=['CHV: \alpha_h=',num2str(alphah,'%3.2f'),' [1/s]   ',...
         '\beta_h=',num2str(betah,'%3.2f'),' [1/s]   ',...
         '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
         'CAV: \alpha_{01}=',num2str(alpha01,'%3.2f'),' [1/s]   ',...
         '\beta_{01}=',num2str(beta01,'%3.2f'),' [1/s]   ',...
         '\beta_{02}=',num2str(beta02,'%3.2f'),' [1/s]   ',...
         '\beta_{03}=',num2str(beta03,'%3.2f'),' [1/s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay,sL,vL,tau,sigma)...
    [x(2);...
     sat(uh(xdelay(1,1),xdelay(2,1),sL(t-tau),vL(t-tau)));...
     x(4);...
     sat(uh(xdelay(3,1),xdelay(4,1),xdelay(1,1),xdelay(2,1)));...
     x(6);...
     sat(u0(xdelay(5,2),xdelay(6,2),xdelay(3,2),xdelay(4,2),xdelay(2,2),vL(t-sigma)))];
 
% leader velocity and position
t0=0; t1=10; t2=30; tend=60;
vL0=vstar; vL1=vstar-10; vL2=vstar; vLend=vstar; 
vL=@(t)(t<t0).*vL0+...
       (t0<=t & t<t1).*(vL0+(vL1-vL0)/(t1-t0)*(t-t0))+...
       (t1<=t & t<t2).*(vL1+(vL2-vL1)/(t2-t1)*(t-t1))+...
       (t2<=t & t<=tend).*(vL2+(vLend-vL2)/(tend-t2)*(t-t2));
sL=@(t)(t<t0).*(vL0*(t-t0))+...
       (t0<=t & t<t1).*(vL0*(t-t0)+(vL1-vL0)/(t1-t0)*(t-t0).^2/2)+...
       (t1<=t & t<t2).*((vL1+vL0)*(t1-t0)/2+...
                        vL1*(t-t1)+(vL2-vL1)/(t2-t1)*(t-t1).^2/2)+...
       (t2<=t & t<=tend).*((vL1+vL0)*(t1-t0)/2+(vL2+vL1)*(t2-t1)/2+...
                        vL2*(t-t2)+(vLend-vL2)/(tend-t2)*(t-t2).^2/2);

% initial condition - history
xinit=@(t)[sL(t)-hstarh-l;vstar;...
           sL(t)-2*hstarh-2*l;vstar;...
           sL(t)-2*hstarh-hstar0-3*l;vstar];

% simulation time
deltat=0.01;
time=t0:deltat:tend;

% perform simulation with built-in DDE23 for solution
sol=dde23(@(t,x,xdelay)model(t,x,xdelay,sL,vL,tau,sigma),[tau,sigma],...
    xinit,[time(1) time(end)]);
x=deval(sol,time);

% extract positions and velocities from the results
sFollow=x(1:2:end,:);
vFollow=x(2:2:end,:);

% calculate headways
hFollow=[sL(time)-sFollow(1,:)-l;sFollow(1:end-1,:)-sFollow(2:end,:)-l];

%% Plot of solution
% plot simulation results for position
figure(1); clf; hold on; box on;
plot(time,sL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,sFollow);
axis([t0,tend,sminplot,smaxplot]);
xlabel('t [s]');
ylabel('s(t) [m]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for headway
figure(2); clf; hold on; box on;
plot([t0,tend],[hstarh,hstarh],'k--');
plot([t0,tend],[hstar0,hstar0],'k--');
set(gca,'ColorOrderIndex',1);
plot(time,hFollow);
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of ',problem,10,parlist]);
legend('h_h*','h_0*','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for speed
figure(3); clf; hold on; box on;
plot([t0,tend],[vstar,vstar],'k--');
plot(time,vL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,vFollow);
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of ',problem,10,parlist]);
legend('v*','CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');