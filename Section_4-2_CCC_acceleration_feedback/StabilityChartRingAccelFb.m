%% Create analytical stability charts for the ring configuration with acceleration feedback
clear; %close all; clc;

%% Parameters
kappa=0.6;
sigma=0.6;
gamma=0.75; % 0, 0.25, 0.5, 0.75
NN=24;
% range of parameters
betamin=-2;         %    0, -2
betamax=3;          %  1.2,  3
alphamin=-1;        %    0, -1
alphamax=4;         %  1.2,  4
aspect_ratio=1;     %    1,  1

% range of frequencies for stability boundaries
ommin=0;
ommax=2*pi;
dom=0.0001;
om=ommin:dom:ommax;

%% Stability boundaries
% omega>0
alpha=zeros(NN-1,length(om));
beta=zeros(size(alpha));
for kk=1:NN-1
    alpha(kk,:)=om.^2.*(cos(om*sigma)-cos(om*sigma-2*kk*pi/NN)+...
                        gamma*(1-cos(2*kk*pi/NN)))./...
                (-om*sin(2*kk*pi/NN)+2*kappa*(1-cos(2*kk*pi/NN)));
    beta(kk,:)=(-om.^2.*(cos(om*sigma)-gamma*cos(2*kk*pi/NN))+...
                    kappa*om.*(sin(om*sigma)-sin(om*sigma-2*kk*pi/NN)-...
                        gamma*sin(2*kk*pi/NN)))./...
                (-om*sin(2*kk*pi/NN)+2*kappa*(1-cos(2*kk*pi/NN)));
    % remove asymptotes
    beta(kk,beta(kk,:)>betamax | beta(kk,:)<betamin)=nan;
end

%% Stability chart
figure(1); clf; hold on; box on;
colors=interp1(linspace(0,1,7),[1,0,0;1,1,0;0,1,0;0,1,1;0,0,1;1,0,1;1,0,0],...
    linspace(0,1,NN));
colormap(colors);
% omega=0
plot([betamin,betamax],[0,0],'Color',colors(1,:));
% omega>0
for kk=1:NN-1
%     plot(beta(kk,:),alpha(kk,:),'b');
    plot(beta(kk,:),alpha(kk,:),'Color',colors(kk+1,:));
end
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the ring configuration with acceleration feedback (N=',num2str(NN),')',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \sigma = ',num2str(sigma,'%3.2f'),' [s]',...
       ' [1/s],   \gamma = ',num2str(gamma,'%3.2f')]);
bar=colorbar; ylabel(bar,'k'); caxis([0,NN]);