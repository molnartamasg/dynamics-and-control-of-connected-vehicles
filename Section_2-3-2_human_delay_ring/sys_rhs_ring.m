%% Right-hand side of equations with N vehicles on a ring
function sys_rhs=sys_rhs_ring(xx,par)
% free parameters
L=par(1);
alpha=par(2);
beta=par(3);
hst=par(5);
vmax=par(6);
kappa=par(7);
amin=par(8);
amax=par(9);
NN=par(10);

% range policy
hgo=hst+vmax/kappa;
V=@(h)vmax*(hgo<=h) + kappa*(h-hst).*(hst<h & h<hgo);
% saturation function
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
% control input
u=@(h0,v0,v1)alpha*(V(h0)-v0)+beta*(v1-v0);

% states are stored as
% xx=[v0(t),v0(t-tau);
%     v1(t),v1(t-tau);
%     ...
%     v{N-1}(t),v{N-1}(t-tau);
%     h0(t),h0(t-tau);
%     h1(t),h1(t-tau);
%     ...
%     h{N-2}(t),h{N-2}(t-tau)];
% note that h{N-1}=L-h0-h1-...-h{N-2}

% right-hand side of equation
v=xx(1:NN,1,:);
vL=[xx(2:NN,1,:);xx(1,1,:)];
vtau=xx(1:NN,2,:);
vLtau=[xx(2:NN,2,:);xx(1,2,:)];
htau=[xx(NN+1:2*NN-1,2,:);L-sum(xx(NN+1:2*NN-1,2,:),1)];
sys_rhs=[sat(u(htau,vtau,vLtau));...
         vL(1:end-1,:,:)-v(1:end-1,:,:)];

end