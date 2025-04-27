%% SetupPOfold - Initialize continuation of folds of periodic orbits
%%
function [pfuncs,pbranch,suc]=SetupPOfold(funcs,branch,ind,varargin)
%% Inputs
%
%  * |funcs|: functions used for DDE
%  * |branch|: branch of psols along which fold was discovered
%  * |ind|: index in points array that is closest to fold (for initial
%      guess)
%
% optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (pbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |pitchfork| (logical, default false) if true then initialization assumes
%   that detected "fold" is pitchfork bifurcation (not safe)
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% All other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
% </html>
%
%% process options
default={'contpar',[],'sys_deri',1e-4,'sys_dtau',1e-4,'correc',true,'dir',[],...
    'step',1e-3,'pitchfork',false,'hjac',1e-4,'df_deriv',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isfield(funcs,'sys_deri_provided') && funcs.sys_deri_provided
    options.sys_deri=funcs.sys_deri;
end
if funcs.tp_del && isfield(funcs,'sys_dtau_provided') && funcs.sys_dtau_provided
    options.sys_dtau=funcs.sys_dtau;
end
branch.point=branch.point(ind); % remove all points but approx fold
%% initialize branch of folds (pbranch) and pass on optional args
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
pfuncs=funcs;
pfuncs.get_comp=@(p,component)extract_from_POfold(p,component,npar);
if ~funcs.tp_del % constant delay
    orig_tau_ind=funcs.sys_tau(); % system delays
    % additional delays needed for extended system:
    [ext_tau,xtau_ind,relations]=generate_tau_ext(point.parameter,orig_tau_ind,...
        'ext_tau_base',period_ind);
    ext_tau_ind=period_ind+(1:length(ext_tau));
    %% set up functions of extended system
    pfuncs.sys_rhs=@(x,p)sys_rhs_POfold(x,p(1:npar),p(beta_ind),p(period_ind),...
        [0,p(orig_tau_ind)],funcs.sys_rhs,dim,xtau_ind,options.sys_deri);
    pfuncs.sys_tau=@()[orig_tau_ind,ext_tau_ind];
    pfuncs.sys_deri=@(x,p,nx,np,v)...
        sys_deri_POfold(x,p,nx,np,v,options.hjac,...
        beta_ind,period_ind,dim,xtau_ind,funcs,pfuncs,options.df_deriv);
    %% required amendments of structures for extended system
    pbranch.parameter.free=[pbranch.parameter.free,beta_ind,period_ind,ext_tau_ind];
else % state-dependent delay
    % indices of additional delays needed for extended system:
    n_orig_tau=funcs.sys_ntau();  % number of state-dependent delays
    xtau_ind=tauSD_ext_ind(n_orig_tau);
    ext_tau=[];
    relations=[];
    %% set up functions of extended system
    pfuncs.sys_rhs=@(x,p)sys_rhs_SD_POfold(x,p(1:npar),p(beta_ind),p(period_ind),...
        funcs.sys_rhs,funcs.sys_tau,options.sys_deri,options.sys_dtau,dim,xtau_ind);
    pfuncs.sys_cond=@(p)sys_cond_POfold(p,funcs.sys_cond,dim,beta_ind,...
        period_ind,pfuncs.get_comp,[]);
    pfuncs.sys_tau=@(itau,x,p)sys_tau_SD_PObif(itau,x,p(1:npar),funcs.sys_tau,dim,xtau_ind);
    pfuncs.sys_ntau=@()n_orig_tau*(n_orig_tau+1);
    pfuncs.sys_deri_provided=false;
    pfuncs.sys_deri=@(x,p,nx,np,v)df_deriv(struct('sys_rhs',pfuncs.sys_rhs),...
        x,p,nx,np,v,options.hjac);
    pfuncs.sys_dtau=@(itau,x,p,nx,np)df_derit(struct('sys_tau',pfuncs.sys_tau),...
        itau,x,p,nx,np,options.hjac);
    %% required amendments of structures for extended system
    pbranch.parameter.free=[pbranch.parameter.free,beta_ind,period_ind];
end
pfuncs.sys_cond=@(p)sys_cond_POfold(p,funcs.sys_cond,dim,beta_ind,...
    period_ind,pfuncs.get_comp,relations);
%% create initial guess for correction
pfoldini0=POfoldInit(funcs,point,branch.method.point,ext_tau,options.pitchfork);
%% correct initial guess and find 2nd point along branch if desired
pbranch.method.point.extra_condition=1;
[pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
    options.dir,options.step,options.correc);
end
