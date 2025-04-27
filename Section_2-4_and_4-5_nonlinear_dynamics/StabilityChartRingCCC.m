%% Create analytical stability charts for the ring setup with time delay
function StabilityChartRingCCC(Ntotal,NN,alphah,betah,kappah,tau,...
                    alpha,kappa,sigma,betamin,betamax,betanmin,betanmax)

% range of frequencies for stability boundaries
ommin=0;
ommax=2*pi;
dom=0.0001;
om=ommin:dom:ommax;

% Human response
n=Ntotal/NN;
Gh=@(s)((betah*s+alphah*kappah)./...
        (s.^2.*exp(s*tau)+(alphah+betah)*s+alphah*kappah)).^(n-1);
Gom=Gh(1i*om); GR=real(Gom); GI=imag(Gom);

% Stability boundaries
% omega>0
betan=zeros(NN-1,length(om));
beta=zeros(size(betan));
for kk=1:NN-1
    D=om.*((GR-1)*sin(2*kk*pi/NN)-GI*(1-cos(2*kk*pi/NN)));
    betan(kk,:)=(-om.^2.*(cos(om*sigma)-GR.*cos(om*sigma-2*kk*pi/NN)-...
        GI.*sin(om*sigma-2*kk*pi/NN))-...
        alpha*om.*(GR*sin(2*kk*pi/NN)+GI*cos(2*kk*pi/NN))+...
        alpha*kappa*(1-2*GR*cos(2*kk*pi/NN)+2*GI*sin(2*kk*pi/NN)+GR.^2+GI.^2))./D;
    beta(kk,:)=(om.^2.*(cos(om*sigma)-cos(om*sigma-2*kk*pi/NN))+...
        alpha*om*sin(2*kk*pi/NN)-...
        alpha*kappa*((GR+1).*(1-cos(2*kk*pi/NN))+GI*sin(2*kk*pi/NN)))./D;
    % remove asymptotes
    beta(kk,beta(kk,:)>betamax | beta(kk,:)<betamin)=nan;
end
% k=0 special case
F=@(om)-om.^2.*((1-real(Gh(1i*om))).*cos(om*sigma)-imag(Gh(1i*om)).*sin(om*sigma))+...
    alpha*kappa*((1-real(Gh(1i*om))).^2+imag(Gh(1i*om)).^2)-...
    alpha*om.*imag(Gh(1i*om));
omstar=arrayfun(@(om0)fsolve(F,om0,optimset('Display','off')),0:pi/10:2*pi);
[~,omidx,~]=unique(round(omstar,4));
omstar=omstar(omidx);
omstar(omstar<1e-2)=[];
betastar=(omstar.^2.*cos(omstar*sigma)-alpha*kappa*(1-real(Gh(1i*omstar))))./...
    (omstar.*imag(Gh(1i*omstar)));

% Stability chart
% figure(1); clf; hold on; box on;
colors=interp1(linspace(0,1,7),[1,0,0;1,1,0;0,1,0;0,1,1;0,0,1;1,0,1;1,0,0],...
    linspace(0,1,NN));
colormap(colors);
% k=0
plot([betastar;betastar],[betanmin;betanmax],'Color',colors(1,:));
% omega>0
for kk=1:NN-1
%     plot(beta(kk,:),betan(kk,:),'b');
    plot(beta(kk,:),betan(kk,:),'Color',colors(kk+1,:));
end
axis([betamin betamax betanmin betanmax]);
% pbaspect([aspect_ratio,1,1]);
xlabel('\beta_{01} [1/s]');
ylabel('\beta_{0n} [1/s]');
title(['Stability chart of the ring setup with time delay',...
       ' (N_{CAV}=',num2str(NN),', N=',num2str(Ntotal),')',10,...
       '\alpha = ',num2str(alpha,'%3.2f'),...
       ' [1/s],   \kappa = ',num2str(kappa,'%3.2f'),...
       ' [1/s],   \sigma = ',num2str(sigma,'%3.2f'),' [s]']);
% bar=colorbar; ylabel(bar,'k'); caxis([0,NN]);

end