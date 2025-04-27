function stability=stst_stabil(funcs,stst,method) 
%% compute spectrum of linearized r.h.s. in equilibrium
% function stability=stst_stabil(stst,method)
% INPUT:
%   funcs problem functions
%	stst steady state point
%	method method parameters 
% OUTPUT:
%	stability stability information
% COMMENT:
%       Assumes (imag(method.lms_parameter_rho)~=0)
%       This condition is tested in the code.
%
% (c) DDE-BIFTOOL v. 3.1.1(27), 14/04/2014
% Added on 05/03/2007 
%
% 
%
%%
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

if imag(method.lms_parameter_rho)==0,
  error(['STST_STABIL: imag(method.lms_parameter_rho)==0 :' ...
	 ' hence, method contains inapproriate paramters for this function. ' ...
	 'Use method=df_mthod(''stst'',1);stst_stabil(stst,' ...
	 ' method.stability), to obtain an appropriate method.stability' ...
	 ' structure for this function.']);
end;
%%
if ~funcs.tp_del
  tau=stst.parameter(sys_tau());
  m=length(tau);
else
  d_ac=method.delay_accuracy;
  m=sys_ntau();
  tau=zeros(1,m);
  xx=repmat(stst.x,[1,m]);
  for j=1:m
    tau(j)=sys_tau(j,xx,stst.parameter);
  end;
  for j=1:m
    if (tau(j)<d_ac),
      s=strcat(['WARNING: delay number_', num2str(j), ...
		' is negative, no stability computed.']);
      disp(s);
      stability=[]; 
      return;
    end;
  end;
end;
taumin=min(tau);
taumax=max(tau);

xx=repmat(stst.x, 1, m+1);
AA=cell(1,m+1);
for j=0:m
  AA{j+1}=sys_deri(xx,stst.parameter,j,[],[]);
end
n=size(AA{1},2);
%%
alpha=method.lms_parameter_alpha;
beta=method.lms_parameter_beta;
k_lms = length(alpha) - 1;

interp_order=method.interpolation_order;
s_min = floor((interp_order-1)/2);
s_plus = (interp_order-1) - s_min;
real_min=method.minimal_real_part;
%% modification by JS to cope with special case of zero delay
if taumax==0 || norm(cat(2,AA{2:end}),'inf')==0
    A=sum(cat(3,AA{:}),3);
    l0=eig(A).';
    if ~isempty(real_min)
        l0=l0(real(l0)>real_min);
    end
    if length(l0)>method.max_number_of_eigenvalues
        l0=l0(1:method.max_number_of_eigenvalues);
    end
    l1=l0.';
    stability=struct('h',NaN,'l0',l0,'l1',l1,'n1',[]);
    return
end
% end of modification (JS)
%%
if isempty(real_min) || (real_min==-Inf) 
  if taumax>0 && taumin>0.1
    real_min=-1/taumin;
  else
    real_min=-1;
  end
end

a_ellipse=real(method.lms_parameter_rho);
b_ellipse=imag(method.lms_parameter_rho);
%% change by JS: adapt nb_nu to keep effort for heuristics constant
% effort ~ nb_nu^ntau/2 * (1+n^3/100) assuming that eig is 100x faster than
% matlab
% approx number of points to create for new heuristics
if isfield(method,'newheuristics_tests')
    npts_n3=method.newheuristics_tests;
else
    npts_n3=2000;
end
nb_nu = floor((npts_n3/(1+n^3/100))^(1/length(tau))); % number of points on unit circle tested
if nb_nu<4 % too few points to approximate circle
    use_newheuristics=false;
else
    use_newheuristics=true;
end
delta_real_min = 0.1; % desired estimated accuracy for roots
%% change by js to avoid error message when heuristics gives no restriction
% (~use_newheuristics | isempty(pts))
h_upperbound=taumax*method.maximal_time_step;
h_lowerbound=taumax*method.minimal_time_step;
if use_newheuristics
    %% point cloud estimating possible spectrum, 
    % original code by K Verheyden, generalized to arbitrary number of
    % delays by JS
    eigA0 = eig(AA{1}).';
    pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0);
else
    pts=[];
end
%% case distinction by JS
if ~isempty(pts)
    %% heuristics gives restriction
    hh=0.9/sqrt(max((real(pts)/a_ellipse).^2+(imag(pts)/b_ellipse).^2));
else
    %% pts can be empty: use old heuristics, according to Verheyden etal
    % J. Comp. Appl. Math. 214 (2007)
    radLMS=min(a_ellipse,b_ellipse); %??
    denominator=abs(real_min)+norm(AA{1},'inf');
    for j=2:m+1
        denominator=denominator+...
            norm(AA{j},'inf')*exp(-abs(real_min)*tau(j-1));
    end
    hh=0.9*radLMS/denominator;
end
%%
hh=min([hh,min(tau(tau>=hh*interp_order*1e-2))/s_plus,h_upperbound]);
hh=max([hh,h_lowerbound]);
% nn=n*(k_lms + ceil(taumax/hh) + s_min);

[mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order);    %#ok<NASGU>

% Note: nn==nL,
% except if some interpolation can be avoided ...

% Throw away mu too close to the origin
ss=exp(real_min*hh);
mu=mu(ss<=abs(mu));

% % Throw away zeros
% mu=find(abs(mu));
% % Zeros were already discarded in the previous step

lambda=mu_to_lambda(mu,hh);

% "A posteriori safeguard", 
% but here we approximate (0.9/h)*LMS^{-1}(T_delta) 
% by (0.9/h)*[ellipse(a_ellipse,b_ellipse)] ...:
lambda=lambda((real(lambda)/a_ellipse).^2+(imag(lambda)/b_ellipse).^2<=(0.9/hh)^2);

[dummy,idx]=sort(real(lambda)); %#ok<ASGLU>
lambda=lambda(idx(end:-1:1));

if length(lambda)>method.max_number_of_eigenvalues
  lambda=lambda(1:method.max_number_of_eigenvalues);
end;

stability=struct('h',hh,'l0',lambda,'l1',[],'n1',[]);

stability=stst_stabil_nwt_corr(stability,AA,tau,method);

return;
