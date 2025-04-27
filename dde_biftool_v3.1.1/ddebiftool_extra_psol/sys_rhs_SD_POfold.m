function y=sys_rhs_SD_POfold(x,p,beta,period,sys_rhs,sys_tau,sys_deri,sys_dtau,dim,xtau_ind)
%% rhs for fold of periodic orbits in SD-DDEs
% xtau_ind(i,:) contains which columns of x(1:dim,:) correspond to
% x(tau(i)+tau(1:end))
% argument sys_deri can be either numeric (say, 1e-3, then a
% finite-difference approximation is used for the linearized system), or a
% function providing the partial derivatives of the DDE's sys_rhs
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
ntau=size(xtau_ind,2)-1;
nvec=size(x,3);
xall=x(1:dim,:,:);
x0=xall(:,xtau_ind(1,:),:);
v=x(dim+1:end,xtau_ind(1,:),:);
y0=sys_rhs(x0,p);
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)sys_rhs(x,p),x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(sys_deri(x0,p,ind-1,[],[]),dev,'*');
end
if isnumeric(sys_dtau) 
    %% no user-provided derivative (sys_dtau is size of deviation)
    dtau=@(x0,dev,itau,ind)app_dir_deriv(@(x)sys_tau(itau,x,p),...
        x0(:,1:itau,:),dev,ind,sys_dtau);
else
    %% sys_dtau is user-provided function
    dtau=@(x0,dev,itau,ind)VAopX(sys_dtau(itau,x0(:,1:itau,:),p,ind-1,[]),dev,'*');
end
%% accumulate xpj=x'(-tau_j), xxj=dxj/dx*v and xTj=dxj/dperiod
on=ones(dim,1);
xx=NaN(dim,ntau+1,nvec);
xT=xx;
xx(:,1,:)=v(:,1,:);
xp=xT;
tau=NaN(1,ntau+1,nvec);
xT(:,1,:)=0;
tau(1,1,:)=0; % count of tau's includes tau0=0
for j=2:ntau+1
    tau(1,j,:)=sys_tau(j-1,x0(:,1:j-1,:),p);
    xp(:,j,:)=sys_rhs(xall(:,xtau_ind(j,:),:),p);
    sumdtau_xx=zeros(1,1,nvec);
    sumdtau_xT=zeros(1,1,nvec);
    for k=1:j-1
        sumdtau_xx=sumdtau_xx+dtau(x0(:,1:j-1,:),xx(:,k,:),j-1,k);
        sumdtau_xT=sumdtau_xT+dtau(x0(:,1:j-1,:),xT(:,k,:),j-1,k);
    end
    xx(:,j,:)=v(:,j,:)-xp(:,j,:).*sumdtau_xx(on,1,:);
    xT(:,j,:)=xp(:,j,:).*(tau(on,j,:)/period-sumdtau_xT(on,1,:));
end

%% accumulate rhs
y1=beta/period*y0;
for j=1:ntau+1
    dev=xx(:,j,:)+beta*xT(:,j,:);
    y1=y1+df(x0(:,1:ntau+1,:),dev,j);
end
y=cat(1,y0,y1);
end
