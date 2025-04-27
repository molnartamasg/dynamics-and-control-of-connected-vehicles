function pfoldini=POfoldInit(funcs,point,method,ext_tau,ispitchfork)
%% crude initial guess for fold of periodic orbits
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
pfoldini=point;
J=psol_jac(funcs,method.collocation_parameters,point.period,...
    point.profile,point.mesh,point.degree,point.parameter,[],1); 
[U,S,V]=svd(J); %#ok<ASGLU>
nullvecs=V(:,end);
v=nullvecs(1:end-1);
if ~ispitchfork
    beta=nullvecs(end);
else
    beta=0;
end
vpoint=p_axpy(0,point,[]);
vpoint.profile=reshape(v,size(point.profile));
vpoint.parameter=beta;
normv=sqrt(p_dot(vpoint,vpoint,'free_par_ind',1));
beta=beta/normv;
vpoint.profile=vpoint.profile/normv;
pfoldini.profile=[pfoldini.profile;vpoint.profile];
pfoldini.parameter=[pfoldini.parameter,beta,pfoldini.period,ext_tau];
end