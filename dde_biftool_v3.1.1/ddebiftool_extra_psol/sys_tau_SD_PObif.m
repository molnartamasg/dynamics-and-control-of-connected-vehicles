function tau=sys_tau_SD_PObif(itau,x,p,sys_tau,dim,xtau_ind)
%% delays of extended systems for fold, torus & period doubling
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
xall=x(1:dim,:,:);
ntau=xtau_ind(1,end)-1;
if itau<=ntau
    tau=sys_tau(itau,xall(:,1:itau,:),p);
else
    [ir,ic]=find(xtau_ind==itau+1,1,'first');
    tau1=sys_tau(ir-1,xall(:,1:xtau_ind(ir,1)-1,:),p);
    tau2=sys_tau(ic-1,xall(:,xtau_ind(ir,1:ic-1),:),p);
    tau=tau1+tau2;
end
end
