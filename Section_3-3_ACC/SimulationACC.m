%% Perform simulations for ACC
clear; %close all; clc;

%% Parameters
% type of vehicle
% vehtype='car';
vehtype='truck';

% vehicle parameters
if strcmp(vehtype,'car')	% personal vehicle
    k=0.45;     % air drag
    m=1500;     % mass
    amin=7;     % max deceleration
    amax=3;     % max acceleration
    P=75000;	% max power
    sigma=0.6;  % delay     % 0.0, 0.6
else	% truck
    k=3.85;     % air drag
    m=29500;	% mass
    amin=6;     % max deceleration
    amax=2;     % max acceleration
    P=300000;	% max power
    sigma=0.8;  % delay
end

% gravitational constant, rolling resistance, elevation
g=9.81;
gamma=0.01;
phi=0.025;

% resistance terms
% p=@(v)0*v;
p=@(v)g*sin(phi)+gamma*g*cos(phi).*sign(v)+k/m*v.^2.*sign(v);
% acceleration constraints
a=@(v)min([amax+0*v;P/m./abs(v)]);
% saturation function
sat=@(v,u)(u<-amin).*(-amin)+(-amin<=u & u<=a(v)).*u+(a(v)<u).*a(v);

% control gains
eta=0.05;	% 0, 0.05
beta=0.5;
alpha=0.4;

% range policy parameters
kappa=0.6;
vmax=30;
hst=5;
hgo=hst+vmax/kappa;

% equilibrium velocity and headway
vstar=20;
hstar=hst+vstar/kappa;

% range policy
V=@(h)vmax*(hgo<=h)+kappa*(h-hst).*(hst<h & h<hgo);
% speed policy
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);

% right-hand side of equations
acc=@(t,x,xsigma,vL,sigma,p,sat,alpha,beta,eta,V,W)...
    [vL(t)-x(2);...
    -p(x(2))+sat(x(2),alpha*(V(xsigma(1))-xsigma(2))+eta*xsigma(3)+beta*(W(vL(t-sigma))-xsigma(2)));...
    V(x(1))-x(2)];

% initial condition - history
xinit=@(t)[hstar;vstar;0];

% leader velocity
vamp=1;
om=2*pi/20;
vL=@(t)vstar+vamp*sin(om*t);

% regions for plotting
hminplot=30;
hmaxplot=50;
vminplot=18;
vmaxplot=22;

% time
t0=0;
tend=100;
deltat=0.01;

%% Simulation
% simulation time
time=t0:deltat:tend;        % simulation time
% perform simulation
if sigma==0
    % use built-in ODE45 for solution
    sol=ode45(@(t,x)acc(t,x,x,vL,sigma,p,sat,alpha,beta,eta,V,W),time,xinit(0));
    x=deval(sol,time);
else
    % use built-in DDE23 for solution
    sol=dde23(@(t,x,xsigma)acc(t,x,xsigma,vL,sigma,p,sat,alpha,beta,eta,V,W),...
                sigma,xinit,[time(1) time(end)]);
    x=deval(sol,time);
end

%% Plot of solution
% plot simulation results for headway
figure(1); clf; hold on; box on;
plot(time,x(1,:),'r');
plot([t0,tend],[hstar,hstar],'k--');
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of ACC with time delay',10,...
        vehtype,'   \sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
        '\eta=',num2str(eta,'%3.2f'),' [1/s^2]   ',...
        '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
        '\alpha=',num2str(alpha,'%3.2f'),' [1/s]']);
legend('h(t)','h*','Location','best');

% plot simulation results for speed
figure(2); clf; hold on; box on;
plot(time,vL(time),'b');
plot(time,x(2,:),'r');
plot([t0,tend],[vstar,vstar],'k--');
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of ACC with time delay',10,...
        vehtype,'   \sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
        '\eta=',num2str(eta,'%3.2f'),' [1/s^2]   ',...
        '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
        '\alpha=',num2str(alpha,'%3.2f'),' [1/s]']);
legend('v_L(t)','v(t)','v*','Location','best');