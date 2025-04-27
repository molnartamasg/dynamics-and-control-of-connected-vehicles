%% plot 2d bifurcation diagram
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
clear
load('LKbifs');
figure(2);clf
ppd=getpar(ind_phi,pdbr);
epd=getpar(ind_eta,pdbr);
ppf=getpar(ind_phi,mwfoldbr);
epf=getpar(ind_eta,mwfoldbr);
plot([-4,8],eta*[1,1],'k-','linewidth',2);
hold on
plot(ph,eh,'ro-',...
    ppd,epd,'ks-','linewidth',2);
plot(ppf,epf,'d-','color',[0,0.5,0],'linewidth',2);
plot( pf,ef,'b.-', pf2,ef2,'b.-', pf2+2*pi,ef2,'b.-','linewidth',2);
axis([-4,8,0,0.009]);
grid on
set(gca,tdeco{:});
legend({'1d bif','RWHopf','MW PD','MW,fold','RWfold'},'location','eastoutside');
xlabel('phi',tdeco{:});
ylabel('eta',tdeco{:});
title('bifs of rotating and modulated waves in Lang-Kobayashi system',tdeco{:});
