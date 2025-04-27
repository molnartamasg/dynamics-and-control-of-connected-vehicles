%% Perform simulations for ACC in mixed traffic
clear; %close all; clc;

%% Parameters
% number of vehicles
N=24;
% distribution of ACC and human driven vehicles
driver=cell(1,N);
n=0; driver(:)={'human'};           % all vehicles are human driven
n=3; driver(n:n:N)={'ACC'};	% replace human driven vehicles with ACC	% n=6, 3, 2

% type of ACC vehicle
vehtype='car';
% vehtype='truck';

% vehicle parameters
if strcmp(vehtype,'car')	% personal vehicle
    k=0.45;     % air drag
    m=1500;     % mass
    amin=7;     % max deceleration
    amax=3;     % max acceleration
    P=75000;	% max power
    sigma=0.6;  % delay
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
phi=0;	% 0.025, 0

% resistance terms
p=@(v)g*sin(phi)+gamma*g*cos(phi).*sign(v)+k/m*v.^2.*sign(v);
% acceleration constraints
a=@(v)min([amax+0*v;P/m./abs(v)]);
% saturation function
sat=@(v,u)(u<-amin).*(-amin)+(-amin<=u & u<=a(v)).*u+(a(v)<u).*a(v);
sath=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;

% ACC parameters
% control gains
eta=0.05;
beta=0.5;
alpha=0.4;

% human driver parameters
alphah=0.1;
betah=0.6;
tau=1.0;

% range policy parameters
kappa=0.6;
vmax=30;
hst=5;
hgo=hst+vmax/kappa;

% time
t0=0;
tend=60;
deltat=0.01;

% leader velocity
vL0=20; vL1=10; vL2=20; 
t1=10; t2=30;
vLead=@(t)(t<t0).*vL0+...
      (t0<=t & t<t1).*(vL0+(vL1-vL0)/(t1-t0)*(t-t0))+...
      (t1<=t & t<t2).*(vL1+(vL2-vL1)/(t2-t1)*(t-t1))+...
      (t2<=t).*vL2;

% range policy
V=@(h)vmax*(hgo<=h)+kappa*(h-hst).*(hst<h & h<hgo);
% speed policy
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);

% right-hand side of equations for ACC vehicle
acc=@(t,x,xdelay,vL,vLdelay,p,sat,alpha,beta,eta,V,W)...
    [vL(t)-x(2);...
    -p(x(2))+sat(x(2),alpha*(V(xdelay(1))-xdelay(2))+eta*xdelay(3)+...
                           beta*(W(vLdelay(t))-xdelay(2)));...
    V(x(1))-x(2)];

% right-hand side of equations for human driven vehicle
human=@(t,x,xdelay,vL,vLdelay,sat,alphah,betah,V)...
    [vL(t)-x(2);...
    sath(alphah*(V(xdelay(1))-xdelay(2))+betah*(vLdelay(t)-xdelay(2)))];

% initial velocity and headway
v0=vLead(t0);
h0=hst+v0/kappa;

% initial condition - history for ACC vehicle
accinit=@(t)[h0;v0;0];
% initial condition - history for human driven vehicle
humaninit=@(t)[h0;v0];

% regions for plotting
vminplot=0;     % 0,  8
vmaxplot=25;    % 25, 22
hminplot=10;	% 10, 20
hmaxplot=55;    % 55, 45

%% Simulation
% simulation time
time=t0:deltat:tend;        % simulation time
% simulate the motion of each vehicle one by one
hFollow=zeros(N,length(time));
vFollow=zeros(size(hFollow));
for kveh=1:N
    disp(['Simulation of vehicle #',num2str(kveh)]);
    % pick initial condition and model for ACC or human driven vehicle
    if strcmp(driver{kveh},'ACC')
        delay=sigma;
        DDErhs=@(t,x,xdelay,vL,vLdelay)acc(t,x,xdelay,vL,vLdelay,p,sat,alpha,beta,eta,V,W);
        DDEic=@(t)accinit(t);
    else
        delay=tau;
        DDErhs=@(t,x,xdelay,vL,vLdelay)human(t,x,xdelay,vL,vLdelay,sat,alphah,betah,V);
        DDEic=@(t)humaninit(t);
    end
    % pick leader velocity for simulations
    if kveh==1  % first vehicle with prescribed velocity
        vL=@(t)vLead(t);
        vLdelay=@(t)vLead(t-delay);
    else        % preceding vehicle
        vL=@(t)interp1(time,vFollow(kveh-1,:),t,'linear','extrap');
        select=@(u,ku)u(ku);
        vLdelay=@(t)(t-delay<=t0).*select(DDEic(t-delay),2)+...
           (t-delay>t0).*interp1(time,vFollow(kveh-1,:),t-delay,'linear','extrap');
    end
    % perform simulation
    sol=dde23(@(t,x,xdelay)DDErhs(t,x,xdelay,vL,vLdelay),delay,DDEic,[time(1) time(end)]);
    x=deval(sol,time);
    % extract headways and velocities from the results
    hFollow(kveh,:)=x(1,:);
    vFollow(kveh,:)=x(2,:);
end

%% Plot of solution
% plot simulation results for headway
figure(3); clf; hold on; box on;
for kv=1:N
    if strcmp(driver(kv),'human')
        plot(time,hFollow(kv,:),'r');
    else
        plot(time,hFollow(kv,:),'b','Linewidth',1.5);
    end
end
plot([t0,tend],[h0,h0],'k--');
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of mixed traffic: human drivers and ACC with delay',10,...
       'ACC: ',vehtype,'   \sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       '\eta=',num2str(eta,'%3.2f'),' [1/s^2]   ',...
       '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]',10,...
       'Human: ','\tau=',num2str(tau,'%3.2f'),' [s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]']);
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:N),driver,'UniformOutput',false);
legend(horzcat(driverlegend,{'h*'}),'Location','northeastoutside');

% plot simulation results for speed
figure(4); clf; hold on; box on;
plot(time,vLead(time),'k');
for kv=1:N
    if strcmp(driver(kv),'human')
        plot(time,vFollow(kv,:),'r');
    else
        plot(time,vFollow(kv,:),'b','Linewidth',1.5);
    end
end
plot([t0,tend],[v0,v0],'k--');
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of mixed traffic: human drivers and ACC with delay',10,...
       'ACC: ',vehtype,'   \sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       '\eta=',num2str(eta,'%3.2f'),' [1/s^2]   ',...
       '\beta=',num2str(beta,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha,'%3.2f'),' [1/s]',10,...
       'Human: ','\tau=',num2str(tau,'%3.2f'),' [s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]']);
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:N),driver,'UniformOutput',false);
legend(horzcat({'#0: leader'},driverlegend,{'v*'}),'Location','northeastoutside');