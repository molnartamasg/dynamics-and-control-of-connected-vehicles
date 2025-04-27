function trini=TorusInit(funcs,point,method,nremove)
%% crude initial guess for stat of toruis bifurcation from Floquet mode
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
[eigval,eigprofile]=mult_crit(funcs,point,method.stability,nremove);
t=repmat(point.mesh,size(point.profile,1),1);
% convert Floquet multiplier mode to Floquet exponent mode (periodic
% function)
eigprofile=eigprofile.*exp(-log(eigval).*t);
omega=atan2(imag(eigval),real(eigval));
upoint=p_axpy(0,point,[]);
upoint.profile=reshape(real(eigprofile),size(point.profile));
vpoint=upoint;
vpoint.profile=reshape(imag(eigprofile),size(point.profile));
utu=p_dot(upoint,upoint);
vtv=p_dot(vpoint,vpoint);
utv=p_dot(upoint,vpoint);
r=1/sqrt(utu+vtv);
gamma=atan2(2*utv,vtv-utu)/2;
qr=r*(upoint.profile*cos(gamma)-vpoint.profile*sin(gamma));
qi=r*(upoint.profile*sin(gamma)+vpoint.profile*cos(gamma));
trini=point;
trini.profile=[trini.profile;qr;qi];
trini.parameter=[trini.parameter,omega/pi,trini.period];
end
