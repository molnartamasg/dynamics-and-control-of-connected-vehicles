%% Create analytical stability charts for the ring setup with time delay
function StabilityChartRing(NN,kappa,tau,betamin,betamax,alphamin,alphamax)

% range of frequencies for stability boundaries
ommin=0;
ommax=2*pi;
dom=0.0001;
om=ommin:dom:ommax;

% Stability boundaries
% omega>0
alpha=zeros(NN-1,length(om));
beta=zeros(size(alpha));
for kk=1:NN-1
    alpha(kk,:)=om.^2.*(cos(om*tau)-cos(om*tau-2*kk*pi/NN))./...
                (-om*sin(2*kk*pi/NN)+2*kappa*(1-cos(2*kk*pi/NN)));
    beta(kk,:)=(-om.^2.*cos(om*tau)+kappa*om.*...
                        (sin(om*tau)-sin(om*tau-2*kk*pi/NN)))./...
                (-om*sin(2*kk*pi/NN)+2*kappa*(1-cos(2*kk*pi/NN)));
    % remove asymptotes
    beta(kk,beta(kk,:)>betamax | beta(kk,:)<betamin)=nan;
end

% Stability chart
% figure(1); clf; hold on; box on;
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
% pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of the ring setup with time delay (N=',num2str(NN),')',10,...
       '\kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \tau = ',num2str(tau,'%3.2f'),' [s]']);
% bar=colorbar; ylabel(bar,'k'); caxis([0,NN]);

end