%% Calculate and plot ReLU and its approximation
clear; clc; % close all;

% constant for ReLU approximation and range for its figure
d=[0.1;1;10];
ymin=-1;
ymax=2;
xmin=-2;
xmax=2;
dx=0.01;
xx=xmin:dx:xmax;

% ReLU function and its approximation
g=@(x)(x+abs(x))/2;
ghat=@(x)(x+x.^2./sqrt(d+x.^2))/2;

% plot ReLU function and its approximation
figure(1); clf; hold on; box on;
plot([0,0],[xmin,xmax],'k--');
plot([xmin,xmax],[0,0],'k--');
plot(xx,g(xx),'k');
plot(xx,ghat(xx));
xlabel('quantity, x'); ylabel('ReLU, g');
axis([xmin,xmax,ymin,ymax]);
pbaspect([(xmax-xmin)/(ymax-ymin),1,1]);