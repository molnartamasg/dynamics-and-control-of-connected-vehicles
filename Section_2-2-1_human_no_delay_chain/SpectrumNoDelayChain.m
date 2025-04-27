%% Create spectrum plots for the vehicle chain
clear; %close all; clc;

%% Parameters
beta=-0.5;	% -0.5, -0.4, 0.5, 0.5,  0.5
alpha=0.4;	%  0.4,  0.4, 0.4,   0, -0.1
kappa=0.6;

% range of plot
Remin=-1;
Remax=1;
Immin=-1;
Immax=1;

%% Eigenvalues
% calculate eigenvalues
s1 = -(alpha+beta)/2 + sqrt((alpha+beta)^2/4 - alpha*kappa);
s2 = -(alpha+beta)/2 - sqrt((alpha+beta)^2/4 - alpha*kappa);

%% Spectrum plot
% plot eigenvalues
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
plot(real(s1),imag(s1),'b*');
plot(real(s2),imag(s2),'b*');
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(s)');
ylabel('Im(s)');
title(['Spectrum of the vehicle chain',10,...
       '\kappa=',num2str(kappa,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f')]);