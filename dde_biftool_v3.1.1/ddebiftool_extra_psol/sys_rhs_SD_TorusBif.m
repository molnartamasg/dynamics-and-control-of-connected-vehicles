function y=sys_rhs_SD_TorusBif(x,p,omega,period,sys_rhs,sys_tau,sys_deri,sys_dtau,dim,xtau_ind)
%% rhs for torus bifurcation of periodic orbits in SD-DDEs
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
ntau=size(xtau_ind,2)-1;
nvec=size(x,3);
xall=x(1:dim,:,:);
x0=xall(:,xtau_ind(1,:),:);
u=x(dim+1:2*dim,xtau_ind(1,:),:);
v=x(2*dim+1:end,xtau_ind(1,:),:);
y0=reshape(sys_rhs(x0,p),[dim,1,nvec]);
rot=pi*omega/period;
if isnumeric(sys_deri) 
    %% no user-provided derivative (sys_deri is size of deviation)
    df=@(x0,dev,ind)app_dir_deriv(@(x)reshape(sys_rhs(x,p),[dim,1,nvec]),...
        x0,dev,ind,sys_deri);
else
    %% sys_deri is user-provided function
    df=@(x0,dev,ind)VAopX(sys_deri(x0,p,ind-1,[],[]),dev,'*');
end
if isnumeric(sys_dtau) 
    %% no user-provided derivative (sys_dtau is size of deviation)
    dtau=@(x0,dev,itau,ind)app_dir_deriv(@(x)reshape(sys_tau(itau,x,p),[1,1,nvec]),...
        x0(:,1:itau,:),dev,ind,sys_dtau);
else
    %% sys_dtau is user-provided function
    dtau=@(x0,dev,itau,ind)VAopX(sys_dtau(itau,x0(:,1:itau,:),p,ind-1,[]),dev,'*');
end
%% accumulate xpj=x'(-tau_j), xxr=dxj/dx*u, xxi=dxj/dx*v
on=ones(dim,1);
xxr=NaN(dim,ntau+1,nvec);
xxi=NaN(dim,ntau+1,nvec);
xxr(:,1,:)=u(:,1,:);
xxi(:,1,:)=v(:,1,:);
xp=xxr;
tau=NaN(1,ntau+1,nvec);
tau(1,1,:)=0; % count of tau's includes tau0=0
for j=2:ntau+1
    tau(1,j,:)=reshape(sys_tau(j-1,x0(:,1:j-1,:),p),[1,1,nvec]);
    xp(:,j,:)=reshape(sys_rhs(xall(:,xtau_ind(j,:),:),p),[dim,1,nvec]);
    sumdtau_xxr=zeros(1,1,nvec);
    sumdtau_xxi=zeros(1,1,nvec);
    for k=1:j-1
        sumdtau_xxr=sumdtau_xxr+dtau(x0(:,1:j-1,:),xxr(:,k,:),j-1,k);
        sumdtau_xxi=sumdtau_xxi+dtau(x0(:,1:j-1,:),xxi(:,k,:),j-1,k);
    end
    c=cos(rot*tau(1,j,:));
    s=sin(rot*tau(1,j,:));
    xxr(:,j,:)= c(on,1,:).*u(:,j,:)+s(on,1,:).*v(:,j,:)-xp(:,j,:).*sumdtau_xxr(on,1,:);
    xxi(:,j,:)=-s(on,1,:).*u(:,j,:)+c(on,1,:).*v(:,j,:)-xp(:,j,:).*sumdtau_xxi(on,1,:);
end

%% accumulate rhs
yu=rot*v(:,1,:);
yv=-rot*u(:,1,:);
for j=1:ntau+1
    udev=xxr(:,j,:);
    yu=yu+df(x0(:,1:ntau+1,:),udev,j);
    vdev=xxi(:,j,:);
    yv=yv+df(x0(:,1:ntau+1,:),vdev,j);
end
y=cat(1,y0,yu,yv);
end
