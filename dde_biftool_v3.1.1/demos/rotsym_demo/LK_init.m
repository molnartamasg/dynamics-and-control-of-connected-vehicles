function [E0,n0,phi0,omega0]=LK_init(alpha,p,eta,tau)
%% obtain initial value for rotating wave for given alpha, pump, eta and tau
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
n0=eta/2;
omega0=n0*alpha-sqrt(eta^2-n0^2);
phi0=angle(n0*(1+1i*alpha)-1i*omega0)+pi+omega0*tau;
phi0=mod(phi0,2*pi);
E0=sqrt((p-n0)/(2*n0+1));
E0=[E0;0];
end
