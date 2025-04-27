%% SetupMWFold - Initialize continuation of folds of modulated waves
%%
function [pfuncs,pbranch,suc]=SetupMWFold(funcs,branch,ind,varargin)
%% Inputs
% 
% * |funcs|: functions used for rotationally symmetric DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind| number of point close to fold
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
% 
% * |contpar| (integers default |[]|): index of continuation parameters 
%   (replacing free parameters in argument branch)
% * |hbif| (default |1e-3|): used for finite differencing when approximating
%   linearized system, replace by |funcs.sys_deri| if |funcs.sys_deri| is
%   analytical
% * |correc| (logical, default |true|): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (|pbranch| has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if dir is non-empty
% * |hjac| (default |1e-6|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%
% (c) DDE-BIFTOOL v. 3.1.1(29), 15/04/2014
%

%% process options
default={'contpar',[],'hbif',1e-3,'correc',true,'dir',[],...
    'step',1e-3,'hjac',1e-6,'df_deriv',true,'minimal_accuracy',1e-5};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind); % remove all points but approx fold
% initialize branch of folds (pbranch)
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters and artificial delays
point=pbranch.point;
if isfield(point,'stability')
    point=rmfield(point,'stability');
end
dim=size(point.profile,1);    % dimension of original problem
npar=length(point.parameter); % number of original system parameters
beta_ind=npar+1;              % location of add. parameter beta
period_ind=npar+2;            % location of copy of period
rho_ind=npar+3;
orig_tau_ind=funcs.sys_tau(); % system delays
% additional delays needed for extended system:
[ext_tau,xtau_ind,relations]=generate_tau_ext(point.parameter,orig_tau_ind,...
    'ext_tau_base',rho_ind);
ext_tau_ind=rho_ind+(1:length(ext_tau));
%% set up functions of extended system
pfuncs=funcs;
pfuncs.get_comp=@(p,component)extract_from_MWFold(p,component,npar);
pfuncs.sys_rhs=@(x,p)sys_rhs_MWFold(x,p(1:npar),p(beta_ind),p(period_ind),...
    p(rho_ind),[0,p(orig_tau_ind)],funcs.sys_rhs,dim,xtau_ind,options.hbif);
pfuncs.sys_cond=@(p)sys_cond_MWFold(p,funcs.sys_cond,...
    dim,beta_ind,period_ind,rho_ind,relations,pfuncs.get_comp);
pfuncs.sys_tau=@()[orig_tau_ind,ext_tau_ind];
pfuncs.sys_deri=@(x,p,nx,np,v)df_deriv(...
    struct('sys_rhs',pfuncs.sys_rhs),x,p,nx,np,v,options.hjac);
%
%pfuncs.sys_deri=@(x,p,nx,np,v)sys_deri_MWFold(x,p,nx,np,v,options.hjac,...
%    beta_ind,period_ind,rho_ind,dim,xtau_ind,funcs,pfuncs,options.df_deriv);
%% required amendments of structures for extended system
pbranch.parameter.free=[pbranch.parameter.free,...
    beta_ind,period_ind,rho_ind,ext_tau_ind];
pbranch.method.point.extra_condition=1;
%% create initial guess for correction
pfoldini0=MWFoldInit(funcs,point,branch.method.point,ext_tau,branch.parameter.free);
%% correct initial guess and find 2nd point along branch if desired
[pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
    options.dir,options.step,options.correc);
end

%% crude initial guess
function pfoldini=MWFoldInit(funcs,point,method,ext_tau,free_par_ind)
pfoldini=point;
J=psol_jac(funcs,method.collocation_parameters,point.period,...
    point.profile,point.mesh,point.degree,point.parameter,free_par_ind,1); 
[rdum,Jcond]=funcs.sys_cond(point); %#ok<ASGLU>
J=[J;Jcond.profile(:)',Jcond.period,Jcond.parameter(free_par_ind)];
[U,S,V]=svd(J); %#ok<ASGLU>
nullvecs=V(:,end);
nx=numel(point.profile);
v=nullvecs(1:nx);
beta=nullvecs(nx+1);
rho=nullvecs(end);
vpoint=p_axpy(0,point,[]);
vpoint.profile=reshape(v,size(point.profile));
vpoint.parameter=[beta,rho];
normv=sqrt(p_dot(vpoint,vpoint,'free_par_ind',[1,2]));
beta=beta/normv;
rho=rho/normv;
vpoint.profile=vpoint.profile/normv;
pfoldini.profile=[pfoldini.profile;vpoint.profile];
pfoldini.parameter=[pfoldini.parameter,beta,pfoldini.period,rho,ext_tau];
end
%% extract components from MWFold
function result_array=extract_from_MWFold(pfold_array,component,npar)
dim=size(pfold_array(1).profile,1)/2;
for i=1:length(pfold_array)
    pfold=pfold_array(i);
    switch component
        case 'kind'
            result='MWFold';
        case 'solution'
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case 'nullvector'
            result=pfold;
            result.profile=result.profile(dim+1:end,:);
            result.parameter=result.parameter(npar+1:3);
        case 'delays'
            result.values=pfold.parameter(npar+4:end);
        case 'omega'
            result=pfold.parameter(npar);
        case 'beta'
            result=pfold.parameter(npar+1);
        case 'rho'
            result=pfold.parameter(npar+3);
        otherwise
            error('MWFold:unknown','component %s unknown',component);
    end
    result_array(i)=result; %#ok<AGROW>
end
end
