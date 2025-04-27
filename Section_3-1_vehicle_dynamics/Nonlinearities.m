%% Calculate and plot nonlinearities related to vehicles
clear; clc; % close all;

% parameters
vmax=30;

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

% plot range
aplotmin=-amin-2;
aplotmax=amax+2;
da=0.01;
aa=aplotmin:da:aplotmax;
vplotmin1=0;
vplotmax1=40;
dv1=0.01;
vv1=vplotmin1:dv1:vplotmax1;
vplotmin2=0;
vplotmax2=60;
dv2=0.01;
vv2=vplotmin2:dv2:vplotmax2;

%% Calculation of the characteristics
% speed policy
W=@(v)(0<v & v<vmax).*v + (vmax<=v)*vmax;

% resistance terms
p=@(v)(g*sin(phi)+gamma*g*cos(phi)+k/m*v.^2).*sign(v);
% acceleration constraints
a=@(v)min([amax+0*v;P/m./v]);

% calculate maximum achievable velocity
vvmax=fsolve(@(v)a(v)-p(v),vmax);
disp(['Maximum velocity: ',num2str(vvmax), ' [m/s]']);

% calculate maximum achievable acceleration
aamax=amax-p(0);
disp(['Maximum achievable acceleartion: ',num2str(aamax), ' [m/s^2]']);

% saturation function
sat=@(u,v)(u<=-amin)*(-amin) + (-amin<u & u<a(v)).*u + (a(v)<=u)*a(v);

%% Plot of the characteristics
% plot speed policy
figure(1); clf; hold on; box on;
plot([vmax,vmax],[vplotmin1,vplotmax1],'k--');
plot([vplotmin1,vplotmax1],[vmax,vmax],'k--');
plot(vv1,W(vv1),'b');
xlabel('speed, v (m/s)'); ylabel('saturated speed, W (m/s)');
axis([vplotmin1,vplotmax1,vplotmin1,vplotmax1]);
pbaspect([1,1,1]);

% plot vehicle characteristics
figure(2); clf; hold on; box on;
plot([vplotmin2,vplotmax2],[0,0],'k--');
LL=zeros(3,1);
LL(1)=plot(vv2,p(vv2),'b');
LL(2)=plot(vv2,a(vv2),'r');
LL(3)=plot(vv2,-amin+0*vv2,'m');
LL(4)=plot(vvmax,a(vvmax),'g.','MarkerSize',24);
xlabel('speed, v [m/s]');
ylabel('acceleration, a [m/s^2]');
axis([vplotmin2,vplotmax2,aplotmin,aplotmax]);
title(['Vehicle characteristics: ',vehtype]);
legend(LL,'p(v)','a_{max}(v)','a_{min}(v)','(v_{max}, a_{max}(v_{max}))',...
    'Location','east');
pbaspect([1,1,1]);

% plot saturation function
figure(3); clf; hold on; box on;
plot([0,0],[aplotmin,aplotmax],'k--');
plot([a(vmax),a(vmax)],[aplotmin,aplotmax],'k--');
plot([-amin,-amin],[aplotmin,aplotmax],'k--');
plot([aplotmin,aplotmax],[0,0],'k--');
plot([aplotmin,aplotmax],[-amin,-amin],'k--');
plot([aplotmin,aplotmax],[a(vmax),a(vmax)],'k--');
plot(aa,sat(aa,vmax),'b');
xlabel('control input, u (m/s^2)'); ylabel('saturated control input, sat (m/s)');
axis([aplotmin,aplotmax,aplotmin,aplotmax]);
pbaspect([1,1,1]);