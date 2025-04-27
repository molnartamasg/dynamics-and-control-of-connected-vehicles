%% Right-hand side of CCC on a ring
function sys_rhs=sys_rhs_ring(xx,par,ind,N,CHVs,CAVs,n,Vh,V0,W,sat)
% free parameters
L=par(ind.L);
alphah=par(ind.alphah);
betah=par(ind.betah);
alpha01=par(ind.alpha01);
beta01=par(ind.beta01);
beta0n=par(ind.beta0n);

% states are stored as
% xx=[v_1(t),v_1(t-sigma),v_1(t-tau);
%     v_2(t),v_2(t-sigma),v_2(t-tau);
%     ...;
%     v_N(t),v_N(t-sigma),v_N(t-tau);
%     h_1(t),h_1(t-sigma),h_1(t-tau);
%     h_2(t),h_2(t-sigma),h_2(t-tau);
%     ...
%     h_{N-1}(t),h_{N-1}(t-sigma),h_{N-1}(t-tau)]

% index of states
hCAVidx = N+CAVs(CAVs<N);
vCAVidx = CAVs;
v1CAVidx = CAVs+1 - (CAVs+1>N)*N;
vnCAVidx = CAVs+n - (CAVs+n>N)*N;
hCHVidx = N+CHVs(CHVs<N);
vCHVidx = CHVs;
v1CHVidx = CHVs+1 - (CHVs+1>N)*N;

% velocities, headways and their derivatives
hCAV = xx(hCAVidx,2,:);
vCAV = xx(vCAVidx,2,:);
v1CAV = xx(v1CAVidx,2,:);
vnCAV = xx(vnCAVidx,2,:);
hCHV = xx(hCHVidx,3,:);
vCHV = xx(vCHVidx,3,:);
v1CHV = xx(v1CHVidx,3,:);
if ismember(N,CAVs)
    hCAV = [hCAV;L - sum(xx([hCAVidx,hCHVidx],2,:),1)];
else
    hCHV = [hCHV;L - sum(xx([hCAVidx,hCHVidx],3,:),1)];
end
hdot = xx([v1CAVidx(CAVs<N),v1CHVidx(CHVs<N)],1,:) - ...
        xx([vCAVidx(CAVs<N),vCHVidx(CHVs<N)],1,:);

% control inputs
uCAV=alpha01*(V0(hCAV)-vCAV)+beta01*(W(v1CAV)-vCAV)+beta0n*(W(vnCAV)-vCAV);
uCHV=alphah*(Vh(hCHV)-vCHV)+betah*(v1CHV-vCHV);

% right-hand side of system
sys_rhs = nan(2*N-1,size(xx,3));
sys_rhs([hCAVidx,hCHVidx],:) = hdot;
sys_rhs(vCAVidx,:) = sat(uCAV);
sys_rhs(vCHVidx,:) = sat(uCHV);
end