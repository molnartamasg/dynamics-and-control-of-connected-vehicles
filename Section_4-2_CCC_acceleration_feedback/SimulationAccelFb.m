%% Perform simulations for CCC with acceleration feedback
clear; %close all; clc;

%% Parameters
% control gains
beta=0.5;
alpha=0.4;
gamma=0.25;	% 0.00, 0.25, 0.50, 0.75
% delay
sigma=0.6;
% acceleration constraints
amin=7;
amax=3;

% range policy parameters
kappa=0.6;
vmax=30;
hst=5;
hgo=hst+vmax/kappa;

% equilibrium velocity and headway
vstar=15;
hstar=hst+vstar/kappa;

% range policy, speed policy, and saturation function
V=@(h)vmax*(hgo<=h)+kappa*(h-hst).*(hst<h & h<hgo);
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);
sat=@(u)max(min(amax,u),-amin);

% control input
u=@(h,v,vL,aL)alpha*(V(h)-v)+beta*(W(vL)-v)+gamma*aL;

% right-hand side of equations
ccc=@(t,x,xdelay,vL,aL,sigma)...
    [vL(t)-x(2);...
    sat(u(xdelay(1),xdelay(2),vL(t-sigma),aL(t-sigma)))];

% initial condition - history
h0=10; v0=vstar;
xinit=@(t)[h0;v0];

% time
t0=0;
tend=8;
deltat=0.01;

% leader velocity
abreak=10;
t1=t0+vstar/abreak;
vL=@(t)vstar*(t<t0)+(vstar-abreak*(t-t0)).*(t0<=t & t<t1);
aL=@(t)-abreak*(t0<=t & t<t1);

% regions for plotting
hminplot=-2;
hmaxplot=12;
vminplot=-2;
vmaxplot=16;
aminplot=-12;
amaxplot=2;

% title to put on figure
problem='CCC with acceleration feedback ';
% list of parameters to put on figure
parlist=['\kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
         '\alpha=',num2str(alpha,'%3.2f'),' [1/s]   ',...
         '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
         '\gamma=',num2str(gamma,'%3.2f'),];

% solver for simulations
solver='ddensd';	% built-in Matlab algorithm
% solver='dde23';     % built-in Matlab algorithm
% solver='ddeab4';	% 4 step Adams-Bashforth method

%% Simulation
% simulation time
time=t0:deltat:tend;        % simulation time
% perform simulation
if strcmp(solver,'ddensd')   % ddensd solver
    % use built-in DDENSD for solution
    sol=ddensd(@(t,x,xdelay,xdotdelay)ccc(t,x,xdelay,vL,aL,sigma),...
                @(t,x)t-sigma,@(t,x)t-sigma,xinit,[time(1) time(end)]);
    [x,xdot]=deval(sol,time);
elseif strcmp(solver,'dde23')   % dde23 solver
    % use built-in DDE23 for solution
    sol=dde23(@(t,x,xdelay)ccc(t,x,xdelay,vL,aL,sigma),...
                sigma,xinit,[time(1) time(end)]);
    [x,xdot]=deval(sol,time);
else	% ab4 solver
    % use own 4 step Adams-Bashforth algorithm
    x=ddeab4(@(t,x,xdelay)ccc(t,x,xdelay,vL,aL,sigma),...
                sigma,xinit,time);
    xdot=gradient(x,1,2)./gradient(time);
end

%% Plot of solution
% plot simulation results for headway
figure(1); clf; hold on; box on;
plot([t0,tend],[hstar,hstar],'k--');
plot([t0,tend],[0,0],'k--');
plot(time,x(1,:),'r');
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of ',problem,10,parlist]);
legend('h*','h=0','h(t)','Location','northeast');

% plot simulation results for speed
figure(2); clf; hold on; box on;
plot([t0,tend],[vstar,vstar],'k--');
plot([t0,tend],[0,0],'k--');
plot(time,vL(time),'b');
plot(time,x(2,:),'r');
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of ',problem,10,parlist]);
legend('v*','v=0','v_L(t)','v(t)','Location','best');

% plot simulation results for acceleration
figure(3); clf; hold on; box on;
plot([t0,tend],[-amin,-amin],'k--');
plot([t0,tend],[amax,amax],'k--');
plot([t0,tend],[0,0],'k--');
plot(time,aL(time),'b');
plot(time,xdot(2,:),'r');
axis([t0,tend,aminplot,amaxplot]);
xlabel('t [s]');
ylabel('a(t) [m/s^2]');
title(['Simulation of ',problem,10,parlist]);
legend('-a_{min}','a_{max}','a=0','a_L(t)','a(t)','Location','best');