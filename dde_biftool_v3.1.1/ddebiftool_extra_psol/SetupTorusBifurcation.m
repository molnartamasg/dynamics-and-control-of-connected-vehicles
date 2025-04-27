%% SetupTorusBifurcation - Initialize continuation of torus or period doubling bifurcations
%%
function [trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which bifurcation was discovered
%
% optional inputs
%
% * |contpar| (integer default |[]|):  set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on
%    torus branch
% * |dir| (integer, default |[]|): which parameter to vary initially along
%    torus branch (trbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to trbranch.method.continuation,
% trbranch.method.point and trbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |trbranch|: bifurcation branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(184), 05/04/2017
% </html>
%
%% process options
default={'contpar',[],'sys_deri',1e-4,'sys_dtau',1e-4,'correc',true,'dir',[],'step',1e-3,...
    'hjac',1e-4,'nremove',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
branch.point=branch.point(ind);
%% initialize branch of torus bifurcations (trbranch) and pass on optional args
trbranch=replace_branch_pars(branch,options.contpar,pass_on);
%% set up numbering and values of additional parameters
point=trbranch.point;
if isfield(point,'stability')
    point=rmfield(point,'stability');
end
dim=size(point.profile,1);    % dimension of original problem
npar=length(point.parameter); % number of original system parameters
omega_ind=npar+1;             % location of add. parameter omega
period_ind=npar+2;            % location of add. parameter (equal to period)
%% set up functions of extended system
trfuncs=funcs;
trfuncs.get_comp=@(p,component)extract_from_tr(p,component);
if ~funcs.tp_del 
    %% constant delays
    tau_ind=funcs.sys_tau();
    trfuncs.sys_rhs=@(x,p)sys_rhs_TorusBif(x,p(1:npar),...
        p(omega_ind),p(period_ind),[0,p(tau_ind)],...
        funcs.sys_rhs,dim,options.sys_deri);
    trfuncs.sys_deri=@(x,p,nx,np,v)...
        sys_deri_TorusBif(x,p,nx,np,v,options.hjac,omega_ind,period_ind,dim,...
        funcs,struct('sys_rhs',trfuncs.sys_rhs));
else
    %% state dependent delay
    n_orig_tau=funcs.sys_ntau();  % number of state-dependent delays
    % additional delays needed for extended system:
    xtau_ind=tauSD_ext_ind(n_orig_tau);
    trfuncs.sys_rhs=@(x,p)sys_rhs_SD_TorusBif(x,p(1:npar),...
        p(omega_ind),p(period_ind),...
        funcs.sys_rhs,funcs.sys_tau,options.sys_deri,options.sys_dtau,dim,xtau_ind);
    trfuncs.sys_ntau=@()n_orig_tau*(n_orig_tau+1);
    trfuncs.sys_deri=@(x,p,nx,np,v)...
        df_deriv(struct('sys_rhs',trfuncs.sys_rhs),x,p,nx,np,v,options.hjac);
    trfuncs.sys_deri_provided=false;
    trfuncs.sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p(1:npar),funcs.sys_tau,dim,xtau_ind);
    %trfuncs.sys_dtau=@(itau,x,p,nx,np)sys_dtau_SD_PObif(itau,x,p,nx,np,funcs.sys_dtau,...
    %    dim,xtau_ind);
    trfuncs.sys_dtau=@(itau,x,p,nx,np)df_derit(struct('sys_tau',trfuncs.sys_tau),...
        itau,x,p,nx,np,options.hjac);
end
trfuncs.sys_cond=@(p)sys_cond_TorusBif(p,funcs.sys_cond,dim,period_ind,...
    trfuncs.get_comp);
%% required amendments of structures for extended system
trbranch.parameter.free=[trbranch.parameter.free,omega_ind,period_ind];
trbranch.method.point.extra_condition=1;
%% create initial guess for correction
trini0=TorusInit(funcs,point,branch.method,options.nremove);
[trbranch,suc]=correct_ini(trfuncs,trbranch,trini0,...
    options.dir,options.step,options.correc);
end
