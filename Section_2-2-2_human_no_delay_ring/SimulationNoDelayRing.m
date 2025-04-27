%% Simulate the car following model

clear; clc; % close all;

%% Parameters
% number of simulated vehicles on ring
N=24;
% driver parameters (identical for each driver)
alpha=0.4;
beta=0.25;  % 0.5, 0.25
hst=5;
hgo=55;
vmax=30;

% simulation time
t0=0; tend=50;
deltat=0.01;
time=t0:deltat:tend;

%% Prescribe initial conditions
% initial velocity of each vehicle
vstar=20;
vamp=5;
kk=1;
om=(alpha+2*beta)*tan(kk*pi/N);
v0=vstar+vamp*cos(2*kk*pi/N*(1:N));
% initial headway of each vehicle
hstar=hst+vstar/vmax*(hgo-hst);     % using equilibrium for linear range policy
h0=hstar+vamp/om*(sin(2*kk*pi/N*(2:N+1))-...
                  sin(2*kk*pi/N*(1:N)));	% perturb headways
% initial state as x0=[h_i(0); v_i(0); h_{i-1}(0); v_{i-1}(0); ...]
x0=[h0;v0]; x0=x0(:);

%% Simulate optimal velocity model of car following
% memory allocation of solution
hFollow=zeros(N,length(time)); % matrix, each row is headway of vehicles
vFollow=zeros(N,length(time)); % matrix, each row is veloctiy of vehicles
% use built-in ODE45 for solution
sol = ode45(@(t,x)OVMNoDelayRing(t,x,alpha,beta,hst,hgo,vmax),...
                        [time(1) time(end)],x0);
xsol = deval(sol,time);
% extract headway and velocity from solution xsol=[h_i; v_i; h_{i-1}; v_{i-1}; ...]
for kveh=1:N
    hFollow(kveh,:)=xsol(2*kveh-1,:).';
    vFollow(kveh,:)=xsol(2*kveh,:).';
end
    
%% Plot simulated data
% plot headway
figure(1); clf; hold on; box on;
plot(time,hFollow);
xlabel('Time [s]');
ylabel('Headway [m]');
title(['Response of ring configuration assuming OVM without delay',10,...
        '\alpha = ',num2str(alpha,'%3.2f'),...
        ' [1/s],   \beta = ',num2str(beta,'%3.2f'),' [1/s]']);
    
% plot velocity
figure(2); clf; hold on; box on;
plot(time,vFollow);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title(['Response of ring configuration assuming OVM without delay',10,...
        '\alpha = ',num2str(alpha,'%3.2f'),...
        ' [1/s],   \beta = ',num2str(beta,'%3.2f'),' [1/s]']);
axis([t0,tend,vstar-2*vamp,vstar+2*vamp]);