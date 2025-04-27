%% Numerical continuation to analyze the nonlinear dynamics of
%% CCC with vehicles on a ring
clear; clc;
% close all;

% add ddebiftool folder to path
addpath(genpath('..\dde_biftool_v3.1.1'));

%% Parameters and system definition
% number of vehicles on ring
N=12;   % 4, 6, 8, 9, 10, 12, 18, 24
% index of CAVs and how many vehicles the CAVs look ahead
n=[]; CAVs = [];      % HV only case
% n=2; CAVs = 1:n:N;	% CHV and CAV case
% n=3; CAVs = 1:n:N;      % CHV and CAV case

% index of CHVs - all other vehicles
CHVs = 1:N; CHVs(CAVs) = [];
% number of CAVs and CHVs
NCAV=length(CAVs);
NCHV=length(CHVs);

% CHV range policy
RPh='cube';	% 'lin', 'quad', 'cube', 'cos'
hsth=5;
hgoh=55;
vmaxh=30;
% CHV parameters
kappah=0.6; % desired kappah - ring length is selected accordingly
if isempty(CAVs)	% HV only case
    alphah=0.4;
    betah=0.5;
    tau=0.6;
else                % CHV and CAV case
    alphah=0.1;
    betah=0.6;
    tau=1;
end

% CAV range policy
RP0='lin';  % 'lin', 'quad', 'cube', 'cos'
hst0=5;
hgo0=55;
vmax0=30;
% CAV parameters
alpha01=0.4;
beta01=0.5;
beta0n=0.1;
sigma=0.6;

% whether to include saturation, and acceleration limits for all vehicles
satur='sat';    % 'sat', 'no_sat'
amin=7;
amax=3;

% range policy of CHV
[Vh,dVh,ddVh]=RP(RPh,vmaxh,hgoh,hsth);
% range policy of CAV
[V0,dV0,ddV0]=RP(RP0,vmax0,hgo0,hst0);

% speed policy of CAV
W=@(vL)vmax0*(vmax0<=vL)+vL.*(vL<vmax0);
dW=@(vL)1*(vL<vmax0);
% W=@(vL)vL;
% dW=@(vL)0*vL+1;

% saturation function for all vehicles
if strcmp(satur,'sat')  % saturation included
    sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
    dsat=@(u)1*(-amin<=u & u<=amax);
else                    % no saturation
    sat=@(u)u;
    dsat=@(u)0*u+1;
end

% ring length based on kappah, with corresponding steady state and kappa0
hstarh=fsolve(@(h)dVh(h)-kappah,hsth/4+hgoh*3/4);
vstar=Vh(hstarh);
hstar0=fsolve(@(h)vstar-V0(h),hstarh);
kappa0=dV0(hstar0);
L=NCHV*hstarh+NCAV*hstar0;

% selected continuation parameters, their range and maximum step size
if isempty(CAVs)    % HV only case
    p1='betah'; p1min=0; p1max=1.2; dp1=0.05; p1label='\beta_{h} [1/s]';
    p2='alphah'; p2min=0; p2max=1.2; dp2=0.05; p2label='\alpha_{h} [1/s]';
else                % CHV and CAV case
    p1='beta01'; p1min=-0.5; p1max=1.5; dp1=0.05; p1label='\beta_{01} [1/s]';
    p2='beta0n'; p2min=-0.5; p2max=1.5; dp2=0.05; p2label='\beta_{0n} [1/s]';
end

% pick a point where periodic solution will be plotted
if isempty(CAVs)    % HV only case
    p1P=0.5;
else                % CHV and CAV case
    p1P=0.25;
end

%% System setup
% vector of initial free parameter values (delays must be included)
parinit=[L,alphah,betah,alpha01,beta01,beta0n,sigma,tau];

% index of free parameters in parameter vector
ind.L=1;        % length of the ring
ind.alphah=2;	% alphah of CHV
ind.betah=3;	% betah of CHV
ind.alpha01=4;	% alpha01 of CAV
ind.beta01=5;	% beta01 of CAV
ind.beta0n=6;	% beta0n of CAV
ind.sigma=7;    % delay of CAV
ind.tau=8;      % delay of HV and CHV

% index of states in state vector
ind.v=1:N;          % velocity
ind.h=N+1:2*N-1;	% headway

% set up right-hand side of system and delays
% funcs=set_funcs('sys_rhs',@(xx,par)sys_rhs_ring(xx,par,ind,N,CHVs,CAVs,n,Vh,V0,W,sat),...
%                 'sys_tau',@()[ind.sigma,ind.tau],'x_vectorized',1);
% set up analytical derivatives of the right-hand side too, if available
funcs=set_funcs('sys_rhs',@(xx,par)sys_rhs_ring(xx,par,ind,N,CHVs,CAVs,n,Vh,V0,W,sat),...
                'sys_deri',@(xx,par,nx,np,v)sys_deri_ring(xx,par,nx,np,v,...
                    ind,N,CHVs,CAVs,n,Vh,dVh,ddVh,V0,dV0,ddV0,W,dW,dsat),...
                'sys_tau',@()[ind.sigma,ind.tau],'x_vectorized',1);

% initial guess for steady state headways and velocities based on range policy
stst0=nan(2*N-1,1);
stst0(1:N)=vstar;               % v*
stst0(N+CAVs(CAVs<N))=hstar0;	% hCAV*
stst0(N+CHVs(CHVs<N))=hstarh;	% hCHV*

%% Branch of equilibria
% set up a branch of steady state solutions
stst_br = SetupStst(funcs,'x',stst0,'parameter',parinit,...
    'contpar',ind.(p1),'max_bound',[ind.(p1),p1max],...
    'min_bound',[ind.(p1),p1min],'max_step',[ind.(p1),dp1/20],'step',dp1/20);

% continue branch
stst_br.method.continuation.plot=0;
stst_br=br_contn(funcs,stst_br,1000);
stst_br=br_rvers(stst_br);
stst_br=br_contn(funcs,stst_br,1000);

% calculate stability along branch
stst_br.method.stability.minimal_real_part=-2;
stst_br=br_stabl(funcs,stst_br,0,0);

% plot stability along branch
figure(2); clf; hold on; box on; grid on;
% plot real part of eigenvalues vs parameter (stability measure)
[xm,ym]=df_measr(1,stst_br);
br_plot(stst_br,xm,ym,'b');
% plot maximum real part vs parameter
ym.subfield='l0';
br_plot(stst_br,xm,ym,'c');
% finalize plot
plot([p1min p1max],[0,0],'r--');
axis([p1min p1max -2 1]);
pbaspect([1 1 1]);
xlabel(p1label); ylabel('\Re(\lambda)');
title('Stability along branch of equilibria');

% locate Hopf bifurcation based on number of unstable exponents (NUE)
stst_br = br_bifdet(funcs, stst_br);  % add Hopf points to the steady state branch
ind_hopf = br_getflags(stst_br,'hopf');
NUE = arrayfun(@(x)sum(real(x.stability.l0)>1e-5),stst_br.point);
% ind_hopf = find(abs(diff(NUE))==2)+1;

%% Branch of periodic solutions
% set up figure for branch
figure(5); clf; hold on; box on; grid on;
xlim([p1min p1max]);
pbaspect([1 1 1]);
xlabel(p1label); ylabel('max(v_0)-min(v_0) [m/s]');
title('Amplitude of periodic solutions');

% resolution of periodic solution over the period
intervals=20;   % 20, 40
degree=3;       %  3,  4
% initial amplitude along branch
radius=0.1;     % 0.01, 0.1

% consider all branches of periodic orbits
psol_branches=cell(length(ind_hopf),1);
for khopf = 1:length(ind_hopf)
    % set up a branch of periodic solutions
    psol_br = SetupPsol(funcs,stst_br,ind_hopf(khopf),...
        'intervals',intervals,'degree',degree,'radius',radius,...
        'contpar',ind.(p1),'max_bound',[ind.(p1),p1max],...
        'min_bound',[ind.(p1),p1min],'max_step',[ind.(p1),dp1/2]);

    % continue branch
    psol_br.method.continuation.plot=1;
    psol_br=br_contn(funcs,psol_br,60);

    % calculate stability along branch
    % psol_br.method.stability.minimal_modulus=0.1;
    psol_br=br_stabl(funcs,psol_br,0,0);

    % store each branch
    psol_branches{khopf}=psol_br;
end

% select branch of periodic orbits - pick a bistable branch if there is any
p1start=cellfun(@(y)y.point(1).parameter(ind.(p1)),psol_branches,...
        'UniformOutput',false);
p1peaks=cellfun(@(y)max(findpeaks(arrayfun(...
        @(x)x.parameter(ind.(p1)),y.point))),psol_branches,...
        'UniformOutput',false);
p1peaks(cellfun(@(y)isempty(y),p1peaks))=p1start(cellfun(@(y)isempty(y),p1peaks));
p1end=cellfun(@(y)y.point(end).parameter(ind.(p1)),psol_branches,...
        'UniformOutput',false);
khopf=find(cell2mat(cellfun(@(x,y,z)x-y>1e-3 & x-z>0,p1peaks,p1start,p1end,...
        'UniformOutput',false)),1,'first');
if isempty(khopf)
    khopf=find(NUE(ind_hopf-1)==0,1,'last');
end
psol_br=psol_branches{khopf};

% plot amplitude, period and stability along branch
figure(6); clf;
p2_st=min(arrayfun(@(x)x.parameter(ind.(p1)),psol_br.point));
p2_e=max(arrayfun(@(x)x.parameter(ind.(p1)),psol_br.point));
% amplitude along branch
subplot(1,3,1); hold on; box on; grid on;
[xm,ym]=df_measr(0,psol_br);
br_plot(psol_br,xm,ym,'b.-');
xlim([p2_st p2_e]);
xlabel(p1label); ylabel('max(v_0)-min(v_0) [m/s]');
title('Amplitude');
% period along branch
subplot(1,3,2); hold on; box on; grid on;
[xm,ym]=df_measr(0,psol_br);
ym.field='period';
ym.col=1;
br_plot(psol_br,xm,ym,'b.-');
xlim([p2_st p2_e]);
xlabel(p1label); ylabel('T_p [s]');
title('Period');
% stability along branch (magnitude of eigenvalues)
subplot(1,3,3); hold on; box on; grid on;
[xm,ym]=df_measr(1,psol_br);
br_plot(psol_br,xm,ym,'b.-');
plot([p2_st,p2_e],[1,1],'r--');
axis([p2_st p2_e 0 1.2]);
xlabel(p1label); ylabel('|\mu|');
title('Eigenvalues');

%% Plot a selected periodic solution
% pick a point and get closest parameter values
p1values=arrayfun(@(x)x.parameter(ind.(p1)),psol_br.point);
% [~,idx]=min(abs(p1values-p1P)); % consider full branch
[~,idx]=min(abs(p1values-p1P)+....
    ([0,diff(p1values)]>=0)*(p1max-p1min));	% consider branch where p1 decreases
psol_p=psol_br.point(idx);

% correct the point to match exactly the selected parameter
psol_p.parameter(ind.(p1))=p1P;
bound_secant=p_axpy(0,psol_br.point(idx),[]);
bound_secant.parameter(ind.(p1))=1;
psol_p=p_correc(funcs,psol_p,ind.(p1),...
    bound_secant,psol_br.method.point,1,psol_br.point(idx));
psol_p.stability=p_stabil(funcs,psol_p,psol_br.method.stability);

% plot the corresponding periodic solution
figure(7); clf;
% velocities
subplot(1,3,1); hold on; box on; grid on;
p_pplot(psol_p,ind.v);
xlabel('t/T_p'); ylabel('v_i [m/s]');
% headways
subplot(1,3,2); hold on; box on; grid on;
p_pplot(psol_p,ind.h);
xlabel('t/T_p'); ylabel('h_i [m]');
title(['Selected periodic solution',...
       ' (',p1,'=',num2str(psol_p.parameter(ind.(p1))),')']);
% eigenvalues (Floquet multipliers)
subplot(1,3,3); hold on; box on; grid on;
p_splot(psol_p);
axis equal;

%% Plot final figure
figure; clf;

% two parameter diagram
subplot(3,2,1); hold on; box on;
if isempty(CAVs)    % HV only case, (betah,alphah) plane
    StabilityChartRing(N,kappah,tau,p1min,p1max,p2min,p2max);
else                % CHV and CAV case, (beta01,beta0n) plane
    StabilityChartRingCCC(N,NCAV,alphah,betah,kappah,tau,...
                    alpha01,kappa0,sigma,p1min,p1max,p2min,p2max);
end
plot([p1min p1max],parinit(ind.(p2))*[1 1],'k');
axis([p1min p1max p2min p2max]);
xlabel(p1label); ylabel(p2label);

% one parameter diagrams
% amplitude
subplot(3,2,3); hold on; box on;
[xm,ym]=df_measr(0,psol_br);
for khopf=1:length(ind_hopf)
    br_plot(psol_branches{khopf},xm,ym,'b-');
end
plot([p1P p1P],[0 30],'k');
plot(p1P,max(psol_p.profile(1,:))-min(psol_p.profile(1,:)),'k*');
axis([p1min p1max 0 30]);
xlabel(p1label); ylabel('max(v_0)-min(v_0) [m/s]');
% period
subplot(3,2,5); hold on; box on;
[xm,ym]=df_measr(0,psol_br);
ym.field='period';
ym.col=1;
for khopf=1:length(ind_hopf)
    br_plot(psol_branches{khopf},xm,ym,'b-');
end
axis([p1min p1max 0 60]);
xlabel(p1label); ylabel('T_p [s]');

% periodic solution
% velocities
subplot(3,2,2); hold on; box on;
plot(psol_p.mesh,psol_p.profile(ind.v,:));
plot([0 1],[vstar vstar],'k--');
axis([0 1 0 35]);
xlabel('t/T_p'); ylabel('v_i [m/s]');
% headways
subplot(3,2,4); hold on; box on;
plot(psol_p.mesh,psol_p.profile(ind.h,:));
plot(psol_p.mesh,psol_p.parameter(ind.L)-sum(psol_p.profile(ind.h,:),1));
if ~isempty(CHVs)
    plot([0 1],[hstarh hstarh],'k--');
end
if ~isempty(CAVs)
    plot([0 1],[hstar0 hstar0],'k--');
end
axis([0 1 0 70]);
xlabel('t/T_p'); ylabel('h_i [m]');