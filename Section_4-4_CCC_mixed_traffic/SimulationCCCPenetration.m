%% Perform simulations for CCC in mixed traffic
clear; %close all; clc;

%% Parameters
% number of vehicles
N=24;
% distribution of CCC and human driven vehicles
driver=cell(1,N);
n=0; driver(:)={'human'};           % all vehicles are human driven
n=3; driver(n:n:N)={'CCC'};	% replace human driven vehicles with CCC	% n=6, 3, 2

% fixed parameters
hst=5;
hgo=55;
vmax=30;
amin=7;
amax=3;
% CHV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;
% CAV
kappa0=0.6;
alpha01=0.4;
beta01=0.5;
beta0n=0.5; % 0.0, 0.5
sigma=0.6;

% range policy
% CHV
Vh=@(h)vmax*(hgo<=h)...
    + kappah*(h-hst).*(hst<h & h<hgo);
% CAV
V0=@(h)vmax*(hgo<=h)...
    + kappa0*(h-hst).*(hst<h & h<hgo);

% saturation function and speed policy
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
W=@(vL)vmax*(vmax<=vL)+vL.*(vL<vmax);

% control input
% CHV
uh=@(h,v,vL)alphah*(Vh(h)-v)+betah*(vL-v);
% CAV
u0=@(h0,v0,v1,vn)alpha01*(V0(h0)-v0)+beta01*(W(v1)-v0)+beta0n*(W(vn)-v0);

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

% right-hand side of equations for CCC vehicle
ccc=@(t,x,xdelay,v1,v1delay,vndelay,sat,u0)...
    [v1(t)-x(2);...
     sat(u0(xdelay(1),xdelay(2),v1delay(t),vndelay(t)))];

% right-hand side of equations for human driven vehicle
human=@(t,x,xdelay,vL,vLdelay,sat,uh)...
    [vL(t)-x(2);...
     sat(uh(xdelay(1),xdelay(2),vLdelay(t)))];

% initial velocity and headway
v0=vLead(t0);
h0h=hst+v0/kappah;
h00=hst+v0/kappa0;

% initial condition - history for CCC vehicle
cccinit=@(t)[h00;v0];
% initial condition - history for human driven vehicle
humaninit=@(t)[h0h;v0];

% regions for plotting
vminplot=8;     % 0,  8
vmaxplot=22;    % 25, 22
hminplot=20;	% 10, 20
hmaxplot=45;    % 55, 45

%% Simulation
% simulation time
time=t0:deltat:tend;        % simulation time
% simulate the motion of each vehicle one by one
hFollow=zeros(N,length(time));
vFollow=zeros(size(hFollow));
for kveh=1:N
    disp(['Simulation of vehicle #',num2str(kveh)]);
    % pick initial condition and model for CCC or human driven vehicle
    if strcmp(driver{kveh},'CCC')
        delay=sigma;
        DDErhs=@(t,x,xdelay,v1,v1delay,vndelay)...
                    ccc(t,x,xdelay,v1,v1delay,vndelay,sat,u0);
        DDEic=@(t)cccinit(t);
    else
        delay=tau;
        DDErhs=@(t,x,xdelay,v1,v1delay,vndelay)...
                    human(t,x,xdelay,v1,v1delay,sat,uh);
        DDEic=@(t)humaninit(t);
    end
    % pick leader velocity for simulations
    select=@(u,ku)u(ku);
    if kveh==1  % first vehicle with prescribed velocity
        v1=@(t)vLead(t);
        v1delay=@(t)vLead(t-delay);
    else        % preceding vehicle
        v1=@(t)interp1(time,vFollow(kveh-1,:),t,'linear','extrap');
        v1delay=@(t)(t-delay<=t0).*select(DDEic(t-delay),2)+...
           (t-delay>t0).*interp1(time,vFollow(kveh-1,:),t-delay,'linear','extrap');
    end
    if kveh<=n
        vndelay=@(t)vLead(t-delay);
    else
        vndelay=@(t)(t-delay<=t0).*select(DDEic(t-delay),2)+...
           (t-delay>t0).*interp1(time,vFollow(kveh-n,:),t-delay,'linear','extrap');
    end

    % perform simulation
    sol=dde23(@(t,x,xdelay)DDErhs(t,x,xdelay,v1,v1delay,vndelay),delay,DDEic,[time(1) time(end)]);
    x=deval(sol,time);
    % extract headways and velocities from the results
    hFollow(kveh,:)=x(1,:);
    vFollow(kveh,:)=x(2,:);
end

%% Plot of solution
% plot simulation results for headway
figure(1); clf; hold on; box on;
for kv=1:N
    if strcmp(driver(kv),'human')
        plot(time,hFollow(kv,:),'r');
    else
        plot(time,hFollow(kv,:),'b','Linewidth',1.5);
    end
end
plot([t0,tend],[h00,h00],'k--');
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of mixed traffic: human drivers and CCC with delay',10,...
       'CCC: ','\sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       '\beta_1=',num2str(beta01,'%3.2f'),' [1/s]   ',...
       '\beta_n=',num2str(beta0n,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha01,'%3.2f'),' [1/s]   ',10,...
       'Human: ','\tau=',num2str(tau,'%3.2f'),' [s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]']);
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:N),driver,'UniformOutput',false);
legend(horzcat(driverlegend,{'h*'}),'Location','northeastoutside');

% plot simulation results for speed
figure(2); clf; hold on; box on;
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
title(['Simulation of mixed traffic: human drivers and CCC with delay',10,...
       'CCC: ','\sigma=',num2str(sigma,'%3.2f'),' [s]   ',...
       '\beta_1=',num2str(beta01,'%3.2f'),' [1/s]   ',...
       '\beta_n=',num2str(beta0n,'%3.2f'),' [1/s]   ',...
       '\alpha=',num2str(alpha01,'%3.2f'),' [1/s]   ',10,...
       'Human: ','\tau=',num2str(tau,'%3.2f'),' [s]   ',...
       '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
       '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]']);
driverlegend=cellfun(@(x,y)['#',num2str(x),': ',y],...
             num2cell(1:N),driver,'UniformOutput',false);
legend(horzcat({'#0: leader'},driverlegend,{'v*'}),'Location','northeastoutside');