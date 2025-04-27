%% Calculate and plot distance-dependent CCC gains
clear; clc; % close all;

% parameters
alpha=0.4;
alphacc=0.9;
beta=0.5;
hgo=30;
hsw=hgo+40;

% plot range
alphaplotmin=0;
alphaplotmax=1;
betaplotmin=0;
betaplotmax=1;
hplotmin=0;
hplotmax=100;
dh=1;
hh=hplotmin:dh:hplotmax;

%% Calculation of the gains
alphabar=@(h)(h<=hsw)*alpha + (h>hsw)*alphacc;
betabar=@(h)(h<=hgo)*beta + (hgo<h & h<=hsw)*beta.*(hsw-h)/(hsw-hgo);

%% Plot of the gains
% alpha gain
figure(1); clf; hold on; box on;
plot([hplotmin,hplotmax],[alpha,alpha],'k--');
plot([hplotmin,hplotmax],[alphacc,alphacc],'k--');
plot([hsw,hsw],[alphaplotmin,alphaplotmax],'k--');
plot(hh,alphabar(hh),'b');
xlabel('headway, h (m)'); ylabel('gain, \alpha (1/s)');
axis([hplotmin,hplotmax,alphaplotmin,alphaplotmax]);
pbaspect([1,1,1]);

% beta gain
figure(2); clf; hold on; box on;
plot([hplotmin,hplotmax],[beta,beta],'k--');
plot([hplotmin,hplotmax],[0,0],'k--');
plot([hgo,hgo],[betaplotmin,betaplotmax],'k--');
plot([hsw,hsw],[betaplotmin,betaplotmax],'k--');
plot(hh,betabar(hh),'b');
xlabel('headway, h (m)'); ylabel('gain, \beta (1/s)');
axis([hplotmin,hplotmax,alphaplotmin,alphaplotmax]);
pbaspect([1,1,1]);