%% Simulate the car following model

clear; clc; % close all;

%% Parameters
% number of simulated follower vehicles
N=24;   % 24, 1
% driver parameters (identical for each driver)
alpha=0.4;
beta=0.5;	% 0.5, 0.25, -0.5
hst=5;
hgo=55;
vmax=30;

%% Prescribe leader velocity profile
% piecewise linear function
vL0=20; vL1=10; vL2=20; vLend=20; 
t0=0; t1=10; t2=30; tend=100;
vL=@(t)(t0<=t & t<t1).*(vL0+(vL1-vL0)/(t1-t0)*(t-t0))+...
       (t1<=t & t<t2).*(vL1+(vL2-vL1)/(t2-t1)*(t-t1))+...
       (t2<=t & t<=tend).*(vL2+(vLend-vL2)/(tend-t2)*(t-t2));
deltah=0;       % perturbation of initial headway
vplotmin=0;     % plot ranges
vplotmax=25;

% % harmonic function
% t0=0; tend=100;
% vL0=20;
% vamp=1;
% omega=0.02*2*pi;
% vL=@(t)vL0+vamp*sin(omega*t);
% deltah=0;       % perturbation of initial headway
% vplotmin=18;    % plot ranges
% vplotmax=22;

% % constant function
% t0=0; tend=20;
% vL0=20;
% vL=@(t)vL0+0*t;
% deltah=1.5;     % perturbation of initial headway
% vplotmin=18;    % plot ranges
% vplotmax=22;

% evaluation of leader's velocity
deltat=0.01;
time=t0:deltat:tend;
vLead=vL(time);

%% Prescribe initial conditions
% initial velocity of each vehicle
vstar=vL0;
v0=vstar*ones(1,N);
% initial headway of each vehicle
hstar=hst+vstar/vmax*(hgo-hst);	% equilibrium for linear range policy
h0=(hstar+deltah)*ones(1,N);
% initial state as x0=[h_i(0); v_i(0); h_{i-1}(0); v_{i-1}(0); ...]
x0=[h0;v0]; x0=x0(:);

%% Simulate optimal velocity model of car following
% memory allocation of solution
hFollow=zeros(N,length(time)); % matrix, each row is headway of vehicles
vFollow=zeros(N,length(time)); % matrix, each row is veloctiy of vehicles
% use built-in ODE45 for solution
sol = ode45(@(t,x)OVMNoDelayChain(t,x,time,vLead,alpha,beta,hst,hgo,vmax),...
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
title(['Response of vehicle chain assuming OVM without delay',10,...
        '\alpha = ',num2str(alpha,'%3.2f'),...
        ' [1/s],   \beta = ',num2str(beta,'%3.2f'),' [1/s]']);
    
% plot velocity
figure(2); clf; hold on; box on;
plot(time,vFollow,time,vLead);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title(['Response of vehicle chain assuming OVM without delay',10,...
        '\alpha = ',num2str(alpha,'%3.2f'),...
        ' [1/s],   \beta = ',num2str(beta,'%3.2f'),' [1/s]']);
axis([t0,tend,vplotmin,vplotmax]);