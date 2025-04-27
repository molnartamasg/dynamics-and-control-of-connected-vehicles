%% Right-hand side of equations with 2 vehicles in an open chain
function sys_rhs=sys_rhs_chain(xx,par)
% free parameters
vstar=par(1);
alpha=par(2);
beta=par(3);
hst=par(5);
vmax=par(6);
kappa=par(7);
amin=par(8);
amax=par(9);

% range policy
hgo=hst+vmax/kappa;
V=@(h)vmax*(hgo<=h) + kappa*(h-hst).*(hst<h & h<hgo);
% saturation function
sat=@(u)(u<-amin).*(-amin)+(-amin<=u & u<=amax).*u+(amax<u).*amax;
% control input
u=@(h0,v0,v1)alpha*(V(h0)-v0)+beta*(v1-v0);

% states are stored as
% xx=[v0(t),v0(t-tau);
%     h0(t),h0(t-tau)];

% right-hand side of equation
sys_rhs=[sat(u(xx(2,2,:),xx(1,2,:),vstar));...
         vstar-xx(1,1,:)];

end