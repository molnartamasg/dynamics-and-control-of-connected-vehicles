function [l,h]=root_app(funcs,x,real_part,alpha,beta,interp_order,par,h_min,h_max,rho,d_ac)
%% approximate rightmost characteristic roots for equilibrium
% function l=root_app(funcs,x,real_part,alpha,beta,interp_order,par,h_min,h_max,rho,d_ac)
% INPUT: 
%   funcs problem functions
%	x steady state in R^n
%	real_part compute roots with real part >= real_part 
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	interp_order order of interpolation in the past
%       par current parameter values in R^p
%	h_min minimal stepsize of discretisation h (relative to max(tau))
%	h_max maximal stepsize of discretisation h (relative to max(tau))
%	rho safety radius of LMS-method
%       d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
% OUTPUT:
%	l approximations of rightmost characteristic roots
%	h stepsize used in discretisation (relative to max(tau))

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

n=length(x);
%% determine delays
if funcs.tp_del==0
  tau=par(sys_tau());
else
  m=sys_ntau();
  xx=x(:,ones(m+1,1));
  tau=NaN(1,m);
  for j=1:m
    tau(j)=sys_tau(j,xx,par);
  end
  for j=1:m,
    if (tau(j)<d_ac),
      s=strcat('WARNING: delay number_',num2str(j),' is negative, no stability computed.');
      disp(s);
      l=[]; 
      h=-1; 
      return
    end
  end
end
%% set boundary for real part depending on max(tau) if not specified
taumax=max(tau);
if isempty(real_part)
  if taumax>0
    real_part=-1/taumax;
  else
    real_part=-1000;
  end
end
%% time step for discrete time map
h=time_h(funcs,x,alpha,beta,real_part,h_min,h_max,par,rho);
%% determine roots
if h>0 % delay case:
  hh=h*taumax;
  S_h=root_int(funcs,x,alpha,beta,hh,interp_order,par);
  nL=size(S_h,2);
  S_h(n+1:nL,1:nL-n)=eye(nL-n);
  s=eig(S_h);
  [dummy,I]=sort(abs(s)); %#ok<ASGLU>
  e=s(I(length(I):-1:1));
  r=exp(real_part*hh);
  l=[];
  for j=1:length(e)
    if abs(e(j))<r
      re=e(1:max(1,j-1));
      l=(log(abs(re))+1i*asin(imag(re)./abs(re)))/hh;
      break
    end
  end
else % ode case: tau(1:m)=0 or D_(1:m)=zeros(n) 
  m=length(tau);
  xx=x(:,ones(m+1,1));
  D=zeros(n);
  for k=0:m
    D=D+sys_deri(xx,par,k,[],[]);
  end;
  le=eig(D);
  l=le(real(le)>=real_part);
end

end
