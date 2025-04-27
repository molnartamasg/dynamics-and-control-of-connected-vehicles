%% Calculate and plot energy contours for CCC in a 4 vehicle scenario
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
sigma=0.6;

% range of parameters
beta02min=-0.5;
beta02max=1.5;
beta03min=-0.5;
beta03max=1.5;
beta02st=40;
beta03st=40;
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
u0=@(s0,v0,s1,v1,v2,v3,beta02,beta03)alpha01*(V0(s1-s0-l)-v0)+...
    beta01*(W(v1)-v0)+beta02*(W(v2)-v0)+beta03*(W(v3)-v0);

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
problem='CCC with 4 vehicles';
% list of parameters to put on figure
parlist0=['CHV: \alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]   ',...
         '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
         '\tau=',num2str(tau,'%3.2f'),' [s]',10,...
         'CAV: \alpha_{01}=',num2str(alpha01,'%3.2f'),' [1/s]   ',...
         '\beta_{01}=',num2str(beta01,'%3.2f'),' [1/s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];

%% Simulation
% right-hand side of equations
model=@(t,x,xdelay,sL,vL,tau,sigma,beta02,beta03)...
    [x(2);...
     sat(uh(xdelay(1,1),xdelay(2,1),sL(t-tau),vL(t-tau)));...
     x(4);...
     sat(uh(xdelay(3,1),xdelay(4,1),xdelay(1,1),xdelay(2,1)));...
     x(6);...
     sat(u0(xdelay(5,2),xdelay(6,2),xdelay(3,2),xdelay(4,2),xdelay(2,2),vL(t-sigma),beta02,beta03))];
 
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
           sL(t)-2*hstarh-2*l;vstar;...
           sL(t)-2*hstarh-hstar0-3*l;vstar];

% simulation time
deltat=0.01;
time=t0:deltat:tend;

% grid in plane of parameters
beta02_v=linspace(beta02min,beta02max,beta02st+1);
beta03_v=linspace(beta03min,beta03max,beta03st+1);
[beta02_m,beta03_m]=ndgrid(beta02_v,beta03_v);

% memory allocation
sFollow=zeros(length(beta02_v),length(beta03_v),length(xinit(0))/2,length(time));
vFollow=zeros(size(sFollow));
hFollow=zeros(size(sFollow));
aFollow=zeros(size(sFollow));
epsFollow=zeros(size(sFollow));
wFollow=zeros(size(sFollow));

% go through all parameter combinations
tic;
for kb02=1:length(beta02_v)
    beta02=beta02_v(kb02);
    for kb03=1:length(beta03_v)
        beta03=beta03_v(kb03);
        % perform simulation with built-in DDE23 for solution
        disp(['Simulation #',num2str((kb02-1)*length(beta03_v)+kb03),'/',...
            num2str(length(beta02_v)*length(beta03_v)),'.']);
        sol=dde23(@(t,x,xdelay)model(t,x,xdelay,sL,vL,tau,sigma,beta02,beta03),...
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
        sFollow(kb02,kb03,:,:)=svalues;
        vFollow(kb02,kb03,:,:)=vvalues;
        hFollow(kb02,kb03,:,:)=hvalues;
        aFollow(kb02,kb03,:,:)=avalues;
        epsFollow(kb02,kb03,:,:)=epsvalues;
        wFollow(kb02,kb03,:,:)=wvalues;
    end
end
toc;

% CAV's final energy consumption
epsCAV=squeeze(epsFollow(:,:,end,end));
wCAV=squeeze(wFollow(:,:,end,end));

% energy optimum
epsmin=min(min(epsCAV));
% wmin=min(min(wCAV(Mmax<1)));
[beta02idx,beta03idx]=find(epsCAV==epsmin);
beta02opt=beta02_v(beta02idx); beta03opt=beta03_v(beta03idx);
wmin=min(min(wCAV));
[beta02lossidx,beta03lossidx]=find(wCAV==wmin);
beta02lossopt=beta02_v(beta02lossidx); beta03lossopt=beta03_v(beta03lossidx);

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
[beta02,beta03,om]=ndgrid(beta02_v,beta03_v,om_v);
[~,~,rho]=ndgrid(beta02_v,beta03_v,rho_v);

% link transfer functions
% HV
T12=(1i*betah*om+alphah*kappah)./...
  (-om.^2.*exp(1i*om*tau)+1i*(alphah+betah).*om+alphah*kappah);
T23=T12;
% CAV
denom=(-om.^2.*exp(1i*om*sigma)+1i*(alpha01+beta01+beta02+beta03).*om+alpha01*kappa0);
T01=(1i*beta01.*om+alpha01*kappa0)./denom;
T02=1i*beta02.*om./denom;
T03=1i*beta03.*om./denom;

% head-to-tail transfer function
G03=T01.*T12.*T23+T02.*T23+T03;
% matrix of maximum head-to-tail transfer function magnitudes
Mmax=max(abs(G03(:,:,0<om_v)),[],3);

% response amplitudes
D=abs(G03).*rho;
% matrix of cost function values
J=sum(D.^2.*om.^2,3);

% cost optimum
Jmin=min(min(J));
[beta02costidx,beta03costidx]=find(J==Jmin);
beta02costopt=beta02_v(beta02costidx); beta03costopt=beta03_v(beta03costidx);

%% Plot of energy contours
% plot energy contours without loss terms
figure(11); clf; hold on; box on;
contour(beta02_m,beta03_m,epsCAV,epsminlevel:2:epsmaxlevel);
% ,'r','ShowText','on'
LL=plot(beta02opt,beta03opt,'kx','Linewidth',1.5);
axis([beta02min,beta02max,beta03min,beta03max]);
colormap(turbo); caxis([epsminlevel,epsmaxlevel]); colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{02} [1/s]');
ylabel('\beta_{03} [1/s]');
title(['Energy contours of ',problem,10,parlist0]);
legend(LL,['\epsilon=',num2str(epsmin),' [J/kg]',10,'(\beta_{02},\beta_{03})=(',...
    num2str(beta02opt),',',num2str(beta03opt),')'],...
    'Location','southwest');

% plot energy contours with loss terms
figure(12); clf; hold on; box on;
contour(beta02_m,beta03_m,wCAV,wminlevel:2:wmaxlevel);
% ,'r','ShowText','on'
LL=plot(beta02lossopt,beta03lossopt,'kx','Linewidth',1.5);
axis([beta02min,beta02max,beta03min,beta03max]);
colormap(turbo); caxis([wminlevel,wmaxlevel]); colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{02} [1/s]');
ylabel('\beta_{03} [1/s]');
title(['Energy contours with loss terms of ',problem,10,parlist0]);
legend(LL,['w=',num2str(wmin),' [J/kg]',10,'(\beta_{02},\beta_{03})=(',...
    num2str(beta02lossopt),',',num2str(beta03lossopt),')'],...
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
% contour(beta02(:,:,1),beta03(:,:,1),Mmax,[1,1],'b','Linewidth',1.5);
contour(beta02(:,:,1),beta03(:,:,1),J,Jminlevel:0.01:Jmaxlevel);
LL=plot(beta02costopt,beta03costopt,'kx','Linewidth',1.5);
axis([beta02min,beta02max,beta03min,beta03max]);
colormap(turbo); caxis([Jminlevel,Jmaxlevel]); c=colorbar;
pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{02} [1/s]');
ylabel('\beta_{03} [1/s]');
title(['Cost contours of ',problem,10,parlist0]);
legend(LL,['J=',num2str(Jmin),' [m^2/s^4]',10,'(\beta_{02},\beta_{03})=(',...
    num2str(beta02costopt),',',num2str(beta03costopt),')'],...
    'Location','southwest');

%% Plot of solution for a selected parameter combination
beta02P=0.5;	% 0.0, 0.5, 0.0, 0.5
beta03P=0.5;	% 0.0, 0.0, 0.5, 0.5
[~,idx1]=min(abs(beta02_v-beta02P));
[~,idx2]=min(abs(beta03_v-beta03P));
parlist=[parlist0,'   \beta_{02}=',num2str(beta02_v(idx1),'%3.2f'),' [1/s]   ',...
                    '\beta_{03}=',num2str(beta03_v(idx2),'%3.2f'),' [1/s]'];

% plot simulation results for position
figure(1); clf; hold on; box on;
plot(time,sL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(sFollow(idx1,idx2,:,:)));
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
plot(time,squeeze(hFollow(idx1,idx2,:,:)));
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
plot(time,squeeze(vFollow(idx1,idx2,:,:)));
axis([t0,tend,vminplot,vmaxplot]);
xlabel('t [s]');
ylabel('v(t) [m/s]');
title(['Simulation of ',problem,10,parlist]);
legend('v*','CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for acceleration
figure(4); clf; hold on; box on;
plot(time,aL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(aFollow(idx1,idx2,:,:)));
axis([t0,tend,aminplot,amaxplot]);
xlabel('t [s]');
ylabel('a(t) [m/s^2]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for energy consumption
figure(5); clf; hold on; box on;
plot(time,epsL(time),'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(epsFollow(idx1,idx2,:,:)));
axis([t0,tend,epsminplot,epsmaxplot]);
xlabel('t [s]');
ylabel('\epsilon(t) [J/kg]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');

% plot simulation results for energy consumption with loss terms
figure(6); clf; hold on; box on;
plot(time,wL,'k');
set(gca,'ColorOrderIndex',1);
plot(time,squeeze(wFollow(idx1,idx2,:,:)));
axis([t0,tend,wminplot,wmaxplot]);
xlabel('t [s]');
ylabel('w(t) [J/kg]');
title(['Simulation of ',problem,10,parlist]);
legend('CHV #3','CHV #2','CHV #1','CAV #0','Location','southeast');

%% Print energy consumption
disp(['Energy consumption with beta02=',num2str(beta02_v(idx1),'%3.2f'),', beta03=',num2str(beta03_v(idx2),'%3.2f'),':']);
disp(['CHV 3:    epsi=',num2str(epsL(time(end)),'%4.1f'),' [J/kg]',9,...
          'Despi=',' 0',' [%]',9,...
          'w=',num2str(wL(end),'%4.1f'),' [J/kg]',9,...
          'Dw=',' 0',' [%]']);
disp(['CHV 2:    epsi=',num2str(epsFollow(idx1,idx2,1,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(idx1,idx2,1,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(idx1,idx2,1,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(idx1,idx2,1,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['CHV 1:    epsi=',num2str(epsFollow(idx1,idx2,2,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(idx1,idx2,2,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(idx1,idx2,2,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(idx1,idx2,2,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['CAV 0:    epsi=',num2str(epsFollow(idx1,idx2,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(idx1,idx2,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(idx1,idx2,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(idx1,idx2,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);

disp('Optima w.r.t. energy (epsilon), energy with loss (w) and cost (J):');
disp(['Oeps:   epsi=',num2str(epsFollow(beta02idx,beta03idx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta02idx,beta03idx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta02idx,beta03idx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta02idx,beta03idx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['Ow:	    epsi=',num2str(epsFollow(beta02lossidx,beta03lossidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta02lossidx,beta03lossidx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta02lossidx,beta03lossidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta02lossidx,beta03lossidx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);
disp(['OJ:     epsi=',num2str(epsFollow(beta02costidx,beta03costidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Despi=',num2str((epsFollow(beta02costidx,beta03costidx,end,end)-epsL(time(end)))/epsL(time(end))*100,'%4.1f'),' [%]',9,...
          'w=',num2str(wFollow(beta02costidx,beta03costidx,end,end),'%4.1f'),' [J/kg]',9,...
          'Dw=',num2str((wFollow(beta02costidx,beta03costidx,end,end)-wL(end))/wL(end)*100,'%4.1f'),' [%]']);