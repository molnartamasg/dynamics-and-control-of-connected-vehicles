%% Calculate and plot energy contours for CCC in a 3 vehicle scenario
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
sigma=0.6;

% range of parameters
beta01min=-0.5;
beta01max=1.5;
beta02min=-0.5;
beta02max=1.5;
beta01st=40;
beta02st=40;
aspect_ratio=1;

% parameters of resistance terms in energy calculation
k=0.45;     % air drag
m=1500;     % mass
gamma=0.01; % rolling resistance
g=9.81;     % gravitational constant
phi=0.0;	% road grade

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

% resistance terms
p=@(v)g*sin(phi)+gamma*g*cos(phi).*sign(v)+k/m*v.^2.*sign(v);

% control input
% CHV
uh=@(s,v,sL,vL)alphah*(Vh(sL-s-l)-v)+betah*(vL-v);
% CAV
u0=@(s0,v0,s1,v1,v2,beta01,beta02)alpha01*(V0(s1-s0-l)-v0)+...
    beta01*(W(v1)-v0)+beta02*(W(v2)-v0);

% regions for plotting
sminplot=-200;
smaxplot=1200;
vminplot=8;
vmaxplot=22;
hminplot=20;
hmaxplot=45;
aminplot=-2;
amaxplot=2;
epsminplot=0;
epsmaxplot=360;
epsminlevel=110;
epsmaxlevel=170;
wminplot=0;
wmaxplot=360;
wminlevel=280;
wmaxlevel=340;
Jminlevel=0.15;
Jmaxlevel=0.75;

% title to put on figure
problem='CCC with 3 vehicles';
% list of parameters to put on figure
parlist0=['CHV: \alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]   ',...
         '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
         '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
         'CAV: \alpha_{01}=',num2str(alpha01,'%3.2f'),' [1/s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay,sL,vL,tau,sigma,beta01,beta02)...
    [x(2);...
     sat(uh(xdelay(1,1),xdelay(2,1),sL(t-tau),vL(t-tau)));...
     x(4);...
     sat(u0(xdelay(3,2),xdelay(4,2),xdelay(1,2),xdelay(2,2),vL(t-sigma),beta01,beta02))];
 
% leader velocity and position
t0=0; t1=10; t2=30; tend=60;
vL0=vstar; vL1=vstar-10; vL2=vstar; vLend=vstar; 
aL0=(vL1-vL0)/(t1-t0); aL1=(vL2-vL1)/(t2-t1); aL2=(vLend-vL2)/(tend-t2);
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
aL=@(t)(t<t0).*0+...
       (t0<=t & t<t1).*aL0+...
       (t1<=t & t<t2).*aL1+...
       (t2<=t & t<=tend).*aL2;
epsL=@(t)(t<t0).*0+...
       (t0<=t & t<t1).*((vL0*(t-t0)+aL0*(t-t0).^2/2).*aL0.*(aL0>0))+...
       (t1<=t & t<t2).*((vL0*(t1-t0)+aL0*(t1-t0).^2/2).*aL0.*(aL0>0)+...
                        (vL1*(t-t1)+aL1*(t-t1).^2/2).*aL1.*(aL1>0))+...
       (t2<=t & t<=tend).*((vL0*(t1-t0)+aL0*(t1-t0).^2/2).*aL0.*(aL0>0)+...
                        (vL1*(t2-t1)+aL1*(t2-t1).^2/2).*aL1.*(aL1>0)+...
                        (vL2*(t-t2)+aL2*(t-t2).^2/2).*aL2.*(aL2>0));

% initial condition - history
xinit=@(t)[sL(t)-hstarh-l;vstar;...
           sL(t)-hstarh-hstar0-2*l;vstar];

% simulation time
deltat=0.01;
time=t0:deltat:tend;

% grid in plane of parameters
beta01_v=linspace(beta01min,beta01max,beta01st+1);
beta02_v=linspace(beta02min,beta02max,beta02st+1);
[beta01_m,beta02_m]=ndgrid(beta01_v,beta02_v);

% memory allocation
sFollow=zeros(length(beta01_v),length(beta02_v),length(xinit(0))/2,length(time));
vFollow=zeros(size(sFollow));
hFollow=zeros(size(sFollow));
aFollow=zeros(size(sFollow));
epsFollow=zeros(size(sFollow));
wFollow=zeros(size(sFollow));

% go through all parameter combinations
tic;
for kb01=1:length(beta01_v)
    beta01=beta01_v(kb01);
    for kb02=1:length(beta02_v)
        beta02=beta02_v(kb02);
        % perform simulation with built-in DDE23 for solution
        disp(['Simulation #',num2str((kb01-1)*length(beta02_v)+kb02),'/',...
            num2str(length(beta01_v)*length(beta02_v)),'.']);
        sol=dde23(@(t,x,xdelay)model(t,x,xdelay,sL,vL,tau,sigma,beta01,beta02),...
            [tau,sigma],xinit,[time(1) time(end)],ddeset('RelTol',1e-6));
        x=deval(sol,time);
        % extract positions and velocities from the results
        svalues=x(1:2:end,:);
        vvalues=x(2:2:end,:);
        % calculate headways
        hvalues=[sL(time)-svalues(1,:)-l;...
                 svalues(1:end-1,:)-svalues(2:end,:)-l];
        % calculate acceleration
        avalues=[vvalues(:,2)-vvalues(:,1),diff(vvalues,1,2)]/deltat;
        % calculation of energy consumption for each point
        epsvalues=cumsum(max(vvalues.*avalues,0),2)*deltat;
        wvalues=cumsum(max(vvalues.*(avalues+p(vvalues)),0),2)*deltat;
        % store results
        sFollow(kb01,kb02,:,:)=svalues;
        vFollow(kb01,kb02,:,:)=vvalues;
        hFollow(kb01,kb02,:,:)=hvalues;
        aFollow(kb01,kb02,:,:)=avalues;
        epsFollow(kb01,kb02,:,:)=epsvalues;
        wFollow(kb01,kb02,:,:)=wvalues;
    end
end
toc;

% CAV's final energy consumption
epsCAV=squeeze(epsFollow(:,:,end,end));
wCAV=squeeze(wFollow(:,:,end,end));

% energy optimum
epsmin=min(min(epsCAV));
% wmin=min(min(wCAV(Mmax<1)));
[beta01idx,beta02idx]=find(epsCAV==epsmin);
beta01opt=beta01_v(beta01idx); beta02opt=beta02_v(beta02idx);
wmin=min(min(wCAV));
[beta01lossidx,beta02lossidx]=find(wCAV==wmin);
beta01lossopt=beta01_v(beta01lossidx); beta02lossopt=beta02_v(beta02lossidx);

% lead vehicle's energy consumption with loss terms
wL=cumsum(vL(time).*(aL(time)+p(vL(time))).*...
            ((vL(time).*(aL(time)+p(vL(time))))>0)*deltat);

%% Cost function based on FFT
% calculate frequency content of leader's velocity
n=2^nextpow2(length(time));
fs=1/deltat;                % sampling frequency
VL=fft(vL(time)-vstar,n)/n; % FFT of leader's velocity fluctuation
VL(1)=VL(1)+vstar;
rho_v=abs(VL(1:n/2+1));     % corresponding amplitudes
rho_v(2:end-1)=2*rho_v(2:end-1);
freq_v=fs*(0:(n/2))/n;      % corresponding frequencies
om_v=2*pi*freq_v;           % corresponding angular frequencies
ommin=min(om_v);
ommax=max(om_v);

% grid in plane of parameters
[beta01,beta02,om]=ndgrid(beta01_v,beta02_v,om_v);
[~,~,rho]=ndgrid(beta01_v,beta02_v,rho_v);

% link transfer functions
% HV
T12=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
% CAV
denom=(-om.^2.*exp(1i*om*sigma)+1i*(alpha01+beta01+beta02).*om+alpha01*kappa0);
T01=(1i*beta01.*om+alpha01*kappa0)./denom;
T02=1i*beta02.*om./denom;

% head-to-tail transfer function
G02=T01.*T12+T02;
% matrix of maximum head-to-tail transfer function magnitudes
Mmax=max(abs(G02(:,:,0<om_v)),[],3);

% response amplitudes
D=abs(G02).*rho;
% matrix of cost function values
J=sum(D.^2.*om.^2,3);

% cost optimum
Jmin=min(min(J));
[beta01costidx,beta02costidx]=find(J==Jmin);
beta01costopt=beta01_v(beta01costidx); beta02costopt=beta02_v(beta02costidx);

%% Plot of energy contours
% plot energy contours without loss terms
figure(11); clf; hold on; box on;
contour(beta01_m,beta02_m,epsCAV,epsminlevel:2:epsmaxlevel);
% ,'r','ShowText','on'
LL=plot(beta01opt,beta02opt,'kx','Linewidth',1.5);
axis([beta01min,beta01max,beta02min,beta02max]);
colormap(turbo); caxis([epsminlevel,epsmaxlevel]); colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{01} [1/s]');
ylabel('\beta_{02} [1/s]');
title(['Energy contours of ',problem,10,parlist0]);
legend(LL,['\epsilon=',num2str(epsmin),' [J/kg]',10,'(\beta_{01},\beta_{02})=(',...
    num2str(beta01opt),',',num2str(beta02opt),')'],...
    'Location','southwest');

% plot energy contours with loss terms
figure(12); clf; hold on; box on;
contour(beta01_m,beta02_m,wCAV,wminlevel:2:wmaxlevel);
% ,'r','ShowText','on'
LL=plot(beta01lossopt,beta02lossopt,'kx','Linewidth',1.5);
axis([beta01min,beta01max,beta02min,beta02max]);
colormap(turbo); caxis([wminlevel,wmaxlevel]); colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{01} [1/s]');
ylabel('\beta_{02} [1/s]');
title(['Energy contours with loss terms of ',problem,10,parlist0]);
legend(LL,['w=',num2str(wmin),' [J/kg]',10,'(\beta_{01},\beta_{02})=(',...
    num2str(beta01lossopt),',',num2str(beta02lossopt),')'],...
    'Location','southwest');

%% Plot of cost function contours
% plot FFT of leader's velocity
figure(21); clf; hold on; box on;
plot(freq_v,rho_v);
axis([0,1,0,vstar]);
xlabel('frequency f [Hz]');
ylabel('amplitude \rho [m/s]');
title('FFT of leader''s velocity');

% plot cost function contours
figure(22); clf; hold on; box on;
% contour(beta01(:,:,1),beta02(:,:,1),Mmax,[1,1],'b','Linewidth',1.5);
contour(beta01(:,:,1),beta02(:,:,1),J,Jminlevel:0.01:Jmaxlevel);
LL=plot(beta01costopt,beta02costopt,'kx','Linewidth',1.5);
axis([beta01min,beta01max,beta02min,beta02max]);
colormap(turbo); caxis([Jminlevel,Jmaxlevel]); c=colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{01} [1/s]');
ylabel('\beta_{02} [1/s]');
title(['Cost contours of ',problem,10,parlist0]);
legend(LL,['J=',num2str(Jmin),' [m^2/s^4]',10,'(\beta_{01},\beta_{02})=(',...
    num2str(beta01costopt),',',num2str(beta02costopt),')'],...
    'Location','southwest');

%% Plot of solution for a selected parameter combination
beta01P=0.5;
beta02P=0.5;    % 0, 0.5
[~,idx1]=min(abs(beta01_v-beta01P));
[~,idx2]=min(abs(beta02_v-beta02P));
parlist=[parlist0,'   \beta_{01}=',num2str(beta01_v(idx1),'%3.2f'),' [1/s]   ',...
                    '\beta_{02}=',num2str(beta02_v(idx2),'%3.2f'),' [1/s]'];

% plot simulation results for position
figure(1); clf; hold on; box on;
plot(time,sL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(sFollow(idx1,idx2,:,:)));
axis([t0,tend,sminplot,smaxplot]);
xlabel('t [s]');
ylabel('s(t) [m]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for headway
figure(2); clf; hold on; box on;
plot([t0,tend],[hstarh,hstarh],'k--');
plot([t0,tend],[hstar0,hstar0],'k--');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(hFollow(idx1,idx2,:,:)));
axis([t0,tend,hminplot,hmaxplot]);
xlabel('t [s]');
ylabel('h(t) [m]');
title(['Simulation of ',problem,10,parlist]);
legend('h_h*','h_0*','CHV #1','CAV #0','Location','southeast');

% plot simulation results for speed
figure(3); clf; hold on; box on;
plot([t0,tend],[vstar,vstar],'k--');
plot(time,vL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(vFollow(idx1,idx2,:,:)));
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of ',problem,10,parlist]);
legend('v*','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for acceleration
figure(4); clf; hold on; box on;
plot(time,aL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(aFollow(idx1,idx2,:,:)));
axis([t0,tend,aminplot,amaxplot]);
xlabel('t [s]');
ylabel('a(t) [m/s^2]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for energy consumption
figure(5); clf; hold on; box on;
plot(time,epsL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(epsFollow(idx1,idx2,:,:)));
axis([t0,tend,epsminplot,epsmaxplot]);
xlabel('t [s]');
ylabel('\epsilon(t) [J/kg]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for energy consumption with loss terms
figure(6); clf; hold on; box on;
plot(time,wL,'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(wFollow(idx1,idx2,:,:)));
axis([t0,tend,wminplot,wmaxplot]);
xlabel('t [s]');
ylabel('w(t) [J/kg]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #2','CHV #1','CAV #0','Location','southeast');

%% Print energy consumption
disp(['Energy consumption with beta01=',num2str(beta01_v(idx1),'%3.2f'),', beta02=',num2str(beta02_v(idx2),'%3.2f'),':']);
disp(['CHV:    epsi=',num2str(epsL(time(end)),'%4.1f'),' [J/kg]',9,...
          'Despi=',' 0',' [%]',9,...
          'w=',num2str(wL(end),'%4.1f'),' [J/kg]',9,...
          'Dw=',' 0',' [%]']);
disp([' HV:    epsi=',num2str(epsFollow(idx1,idx2,1,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(idx1,idx2,1,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(idx1,idx2,1,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(idx1,idx2,1,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['CAV:    epsi=',num2str(epsFollow(idx1,idx2,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(idx1,idx2,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(idx1,idx2,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(idx1,idx2,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);

disp('Optima w.r.t. energy (epsilon), energy with loss (w) and cost (J):');
disp(['Oeps:   epsi=',num2str(epsFollow(beta01idx,beta02idx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta01idx,beta02idx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta01idx,beta02idx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta01idx,beta02idx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['Ow:	    epsi=',num2str(epsFollow(beta01lossidx,beta02lossidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta01lossidx,beta02lossidx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta01lossidx,beta02lossidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta01lossidx,beta02lossidx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['OJ:     epsi=',num2str(epsFollow(beta01costidx,beta02costidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta01costidx,beta02costidx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta01costidx,beta02costidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta01costidx,beta02costidx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);