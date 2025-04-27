%% Create spectrum plots for the vehicle chain with ZOH
clear; %close all; clc;

%% Parameters and system definition
% following vehicle
kappa=0.6;
deltat=0.4;
beta=-0.5;	% -0.5,  0.5,  0.5
alpha=0.4;	%  0.4,  0.4, -0.1

% range of plot
Remin=-1.2;
Remax=1.2;
Immin=-1.2;
Immax=1.2;

%% First equilibrium point
% expressions of the transfer function
syms k dt a b z om real
T=((dt^2/2*a*k+dt*b)*z+dt^2/2*a*k-dt*b)/...
    (z^3-2*z^2+(dt^2/2*a*k+dt*(a+b)+1)*z+dt^2/2*a*k-dt*(a+b));
Tom=subs(T,z,exp(1i*om*dt));
[N,D]=numden(T);

% get corresponding eigenvalues
charpol=subs(D,[k,dt,a,b],[kappa,deltat,alpha,beta]);
eig_v=double(root(charpol));

% plot its eigenvalues
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
plot(cos(0:2*pi/100:2*pi),sin(0:2*pi/100:2*pi),'k--');	% unit circle
plot(real(eig_v),imag(eig_v),'r*','Linewidth',1);
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(z)');
ylabel('Im(z)');
title(['Spectrum of the vehicle chain with ZOH',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \Deltat=',num2str(deltat,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);