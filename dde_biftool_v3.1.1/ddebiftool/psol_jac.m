function [J,res,tT,extmesh]=psol_jac(funcs,c,T,psol_prof,t,deg,par,free_par,ph,varargin)
%% residual & Jacobian of collocation problem for periodic orbits
% function [J,res,tT,extmesh]=psol_jac(c,T,profile,t,deg,par,free_par,phase,varargin)
% INPUT:
%   funcs problem functions
%	c collocation parameters in [0,1]^m
%	T period
%	profile profile in R^(n x deg*l+1)
%	t representation points in [0,1]^(deg*l+1)
%	deg degree piecewise polynomial
%	par current parameter values in R^p
%	free_par free parameters numbers in N^np
%	phase use phase condition or not (s = 1 or 0)
%   wrapJ (optional key-value pair, default true) 
%        wrap time points periodically into [0,1]
% OUTPUT:
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+np)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% modified to permit arbitrary nesting of state-dependent delays,
% vectorisation and optional re-use for computation of Floquet multipliers
% and extra conditions for state-dependent delay (p_tau, etc)
%
% Optional inputs:
%
% If 'wrapJ' (default true) the returned jacobian is augmented with
% derivative wrt period and free_par and wrapped around periodically inside
% [0,1] if ~wrapJ the returned jacobian puts its entries into the interval
%  [-min(delay,0),max(1-[0,delays])], no augmentation is done. This
%  Jacobian can be used for computation of Floquet multipliers and modes
%
% If 'c_is_tvals' (default false) then the values in c (must be not empty) are
% taken as those values of time in [0,1] of the entire period, where
% residual and Jacobian are calculated.
%
% The argument 'Dtmat' (default=Id) is multiplied with the time
% derivative. Dtmat==zeros() permits algebraic equations.
%
% The argument 'bc' (default true) controls whether boundary conditions
% should be attached.

% (c) DDE-BIFTOOL v. 3.1.1(106), 22/08/2015
%
% 
%
%% optional
default={'wrapJ',true,'c_is_tvals',false,'bc',true,'Dtmat',eye(size(psol_prof,1)),...
    'rotationcheck',true};
options=dde_set_options(default,varargin,'pass_on');
%% use functions from funcs
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;

%% define problem & problem dimensions
n=size(psol_prof,1);      % dimension of x
nf=size(options.Dtmat,1); % dimension of f
np=length(free_par);      % number of free parameters
tp_del=funcs.tp_del;
if tp_del==0         % constant delay
    n_tau=sys_tau(); % delay numbers
    tau=par(n_tau);  % delay values
    tau=[0;tau(:)];  % incl tau=0
    d=length(n_tau)+1; % number of delays
    tau_dep=[];
else                 % state dependent delay
    if isfield(funcs,'ntau_is_matrix')
        %% check if dependencies of tau(xx(:,k) are explicitly provided
        tau_dep=sys_ntau();
        d=size(tau_dep,1)+1;
        tau_dep=[zeros(1,size(tau_dep,2));tau_dep];
    else
        %% assume tau(i) depends on xx(:,1:i)
        d=sys_ntau()+1; % number of delays
        tau_dep=false(d);
        for i=1:d-1
            tau_dep(i+1,1:i)=true;
        end
    end
end

nint=(length(t)-1)/deg; % number of collocation intervals
if ~options.c_is_tvals
    neqs=deg*nint;          % number of equations from collocation points
else
    neqs=length(c); % only evaluate res,Jac at times c
    ph=false;
    options.bc=false;
end
%% check array sizes:
if nint~=floor(nint)
    error('PSOL_JAC: t (length%d) does not contain %d intervals of %d points!',...
        length(t),nint,deg);
end
if ~options.c_is_tvals && length(c)~=deg && ~isempty(c)
    error('PSOL_JAC: wrong number of collocation parameters length(c)=%d, m=%d!',...
        length(c),deg);
end
tcoarse=t(1:deg:end); % boundaries of collocation intervals
h_int=diff(tcoarse); % lengths of collocation intervals
%% obtain collocation parameters and times
if ~options.c_is_tvals
    if isempty(c)
        c=poly_gau(deg);
        c=c(:)';
    else
        c=c(:)';
    end
    t_c=tcoarse(ones(length(c),1),1:end-1)+...
        c(ones(length(h_int),1),:)'.*h_int(ones(length(c),1),:);
else
    t_c=c;
end
t_c=t_c(:)'; % evaluation points for DE residuals
%% pre-allocate needed temporary arrays
Pb=zeros(d,deg+1,neqs);  % Lagrange stamp, mapping delayed values onto profile
dPb=zeros(d,deg+1,neqs); % Lagrange stamp, mapping delayed derivatives onto profile
xx=zeros(n,d,neqs);  % array [x(t),x(t-tau1),...]
dtx=zeros(n,d,neqs);  % array [x'(t),x'(t-tau1),...]
c_tau=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t
c_tau_mod=zeros(d,neqs); % time points t-tau(i-1) for current colloc point t (mod[0,1])
index_b=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1)
index_b_mod=zeros(d,neqs); % starting indices of collocation intervals for t-tau(i-1) (mod[0,1])
hhh=zeros(d,neqs); % length of collocation interval in which t-tau(i-1) lies
c_tau_trans=zeros(d,neqs); % position of t-tau(i) in colloc interval, rescaled to [0,1]
if tp_del~=0
    tT=NaN(d,neqs);
    tT(1,:)=0;
else
    tT=tau(:)/T;
    tT=repmat(tT,1,neqs);
end
%% for all collocation points find delays, delayed profiles and derivatives
% (vectorization in l_i & m_i possible)
for t_i=1:d
    c_tau(t_i,:)=t_c-tT(t_i,:);
    c_tau_mod(t_i,:)=mod(c_tau(t_i,:),1);
    [xdum,ibase,c_tau_trans(t_i,:)]=psol_eva(t,t,c_tau_mod(t_i,:),deg); %#ok<ASGLU>
    index_b_mod(t_i,:)=(ibase-1)*deg+1;
    index_b(t_i,:)=index_b_mod(t_i,:)+round(c_tau(t_i,:)-c_tau_mod(t_i,:))*nint*deg;
    hhh(t_i,:)=h_int(ibase);
    Pb(t_i,:,:)=poly_elg(deg,c_tau_trans(t_i,:));
    dPb(t_i,:,:)=reshape(poly_del(deg,c_tau_trans(t_i,:)),[deg+1,neqs])./...
        hhh(t_i*ones(deg+1,1),:);
    for k=1:neqs
        xx(:,t_i,k)=psol_prof(:,index_b_mod(t_i,k)+(0:deg))*Pb(t_i,:,k)';
        dtx(:,t_i,k)=psol_prof(:,index_b_mod(t_i,k)+(0:deg))*dPb(t_i,:,k)';
        if tp_del~=0 && t_i<d
            tT(t_i+1,k)=sys_tau(t_i,xx(:,1:t_i,k),par)/T;
        end
    end
end
%% determine index ranges for wrapped or unwrapped Jacobian
if ~options.wrapJ % compute Jde_dx only for Floquet multipliers
    indmin=min(index_b(:));
    indshift=1-indmin;
    indmax=max(index_b(:))+deg;
    doJcomb=false;
    % set up extended mesh
    indrg=indmin:indmax;
    t_ind=mod(indrg-1,nint*deg)+1;
    t_shift=floor((indrg-1)/(nint*deg));
    extmesh=t(t_ind)+t_shift;
    index_b=index_b+indshift;
else % compute augmented jacobian
    indshift=0;
    indmax=deg*nint+1;
    index_b=index_b_mod;
    doJcomb=true;
    extmesh=t;
end
%% init J, res:
resde=zeros(nf,neqs);
Jde_dx=zeros(nf,neqs,n,indshift+indmax);
if doJcomb
    Jde_dT=zeros(nf,neqs);
    Jde_dp=zeros(nf,neqs,np);
    if options.bc
        Jbc_dx=zeros(n,n,neqs+1);
        Jbc_dp=zeros(n,np);
        Jbc_dT=zeros(n,1);
    end
end
%% obtain all values of sys_rhs, sys_deriv and sys_dxtau
vals=psol_sysvals(funcs,xx,par,free_par,tp_del,tau_dep,'fdim',nf);
%%
for i=1:neqs
    %% precompute all Jacobians for delayed values after delayed values are known
    if tp_del~=0
        dxxdx=zeros([n,n,d,d]); %% dxx(:,k)/dx (sdddes)
        dxxdp=zeros([n,np,d]); %% dxx(:,k)/dp (sdddes)
        dxxdT=zeros([n,d]); %% dxx(:,k)/dT (sdddes)
        dxxdx(:,:,1,1)=eye(n);
        for t_i=2:d
            sel=find(tau_dep(t_i,:));
            dxxdx(:,:,t_i,t_i)=eye(n);
            prefac=dtx(:,t_i,i)/T;
            dxxdp(:,:,t_i)=-prefac*vals.dtau_dp(t_i,:,i);
            dxxdT(:,t_i)=prefac*tT(t_i,i);
            for t_k=sel
                fac=prefac*vals.dtau_dx(:,t_k,t_i,i)';
                dxxdxk=reshape(fac*reshape(dxxdx(:,:,:,t_k),[n,n*d]),[n,n,d]);
                dxxdx(:,:,:,t_i)=dxxdx(:,:,:,t_i)-dxxdxk;
                dxxdp(:,:,t_i)=dxxdp(:,:,t_i)-fac*dxxdp(:,:,t_k);
                dxxdT(:,t_i)=dxxdT(:,t_i)-fac*dxxdT(:,t_k);
            end
        end
    end
    %% insert all derivatives into large Jacobians
    %% add dtx for x' in Jx and res:
    resde(:,i)=options.Dtmat*dtx(:,1,i);
    Jde_dx(:,i,:,index_b(1,i)+(0:deg))=Jde_dx(:,i,:,index_b(1,i)+(0:deg))...
        +reshape(kron(dPb(1,:,i),options.Dtmat),[nf,1,n,deg+1]);
    if doJcomb
        %% add parameter in Jde_dp:
        for p_i=1:np
            Jde_dp(:,i,p_i)=Jde_dp(:,i,p_i)-T*vals.dfdp(:,p_i,i);
        end
        %% add -f for dT in J:
        Jde_dT(:,i)=Jde_dT(:,i)-vals.f(:,i);
    end
    %% add Tf in res:
    resde(:,i)=resde(:,i)-T*vals.f(:,i);
    %% delayed values (incl tau=0): insert Jacobians into J
    for t_i=1:d
        if tp_del==0
            %% add -T*Pb*dfdx(:,t_i)
            Jde_dx(:,i,:,index_b(t_i,i)+(0:deg))=...
                Jde_dx(:,i,:,index_b(t_i,i)+(0:deg))-...
                reshape(kron(Pb(t_i,:,i),T*vals.dfdx(:,:,t_i,i)),[nf,1,n,deg+1]);
            if doJcomb
                %% add -T*A1*sum b*dP*dc_tau for dT in J:
                Jde_dT(:,i)=Jde_dT(:,i)-vals.dfdx(:,:,t_i,i)*dtx(:,t_i,i)*tT(t_i,i);
            end
        else
            for t_k=1:t_i
                dfdxxik=T*vals.dfdx(:,:,t_i,i)*dxxdx(:,:,t_k,t_i);
                Jde_dx(:,i,:,index_b(t_k,i)+(0:deg))=...
                    Jde_dx(:,i,:,index_b(t_k,i)+(0:deg))-...
                    reshape(kron(Pb(t_k,:,i),dfdxxik),[nf,1,n,deg+1]);
            end
            if doJcomb
                Jde_dp(:,i,:)=Jde_dp(:,i,:)-reshape(T*vals.dfdx(:,:,t_i,i)*dxxdp(:,:,t_i),[nf,1,np]);
                Jde_dT(:,i)=Jde_dT(:,i)-T*vals.dfdx(:,:,t_i,i)*dxxdT(:,t_i);
            end
        end
        %% parameter derivative if delay is free par and ~sd-dde
        if tp_del==0 && doJcomb && t_i>1
            p_i=find(free_par==n_tau(t_i-1),1,'first');
            if ~isempty(p_i)
                Jde_dp(:,i,p_i)=Jde_dp(:,i,p_i)+vals.dfdx(:,:,t_i,i)*dtx(:,t_i,i);
            end
        end
    end
end
%% assemble overall residual & Jacobian
Jde_dx=reshape(Jde_dx,[nf*neqs,n*(indshift+indmax)]);
res=resde(:);
J=Jde_dx;
if ~doJcomb
    return
end
J=[J,Jde_dT(:),reshape(Jde_dp,[nf*neqs,np])];
%% periodicity condition:
if options.bc
    Jbc_dx(:,:,1)=eye(n);
    Jbc_dx(:,:,end)=-eye(n);
    resbc=psol_prof(:,1)-psol_prof(:,size(psol_prof,2));
    if options.rotationcheck
        resbc=mod(resbc+pi,2*pi)-pi;
    end
    res=[res;resbc(:)];
    J=[J;reshape(Jbc_dx,[n,n*(neqs+1)]),Jbc_dT,Jbc_dp];
end
%% phase condition:
if ph 
    point=struct('kind','psol','parameter',par,'mesh',t,'degree',deg,...
        'profile',psol_prof,'period',T);
    p0=p_axpy(0,point,[]);
    [resph,p_Jph_dx]=p_dot(p0,point,'derivatives',[0,1]);
    Jph_dx=p_Jph_dx.profile(:)';
    res=[res;resph];
    J=[J; Jph_dx(:)',0,zeros(1,np)];
end
end
