function S=root_int(funcs,x,alpha,beta,h,interp_order,par)
%% first n rows of integration operator S(h) in R^(n x L)
% function S=root_int(funcs,x,alpha,beta,h,interp_order,par)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	h absolute stepsize > 0
%	interp_order order of interpolation in the past
%       par current parameter values in R^p
% 
% OUTPUT:
%	S first n rows of integration operator S(h) in R^(n x L)

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

if funcs.tp_del==0
  tau=par(sys_tau());
  m=length(tau);
  xx=x(:,ones(m+1,1));
else
  m=sys_ntau();
  xx=x;
  tau=NaN(1,m);
  for j=1:m
    tau(j)=sys_tau(j,xx,par);
  end;
end;

n=length(x);
k=length(alpha);

% present terms:

A=sys_deri(xx,par,0,[],[]);

if abs(beta(k))>0
  fac=inv(alpha(k)*eye(n)-h*beta(k)*A);
else
  fac=eye(n)/alpha(k);
end;
S=zeros(n,n*(k-1));
for j=1:k-1
  S(1:n,(j-1)*n+(1:n))=fac*(-alpha(k-j)*eye(n)+h*beta(k-j)*A);
end;

% past terms:

interp_degree=interp_order-1;

s=ceil(interp_degree/2);
r=floor(interp_degree/2);

for p=1:m
  B=sys_deri(xx,par,p,[],[]);
  l=ceil(tau(p)/h);
  if s>l-2
    ss=l-2;
    rr=interp_order-ss;
  else
    ss=s;
    rr=r;
  end;
  epsi=l-tau(p)/h;
  gamma_vect=time_nrd(epsi,rr,ss);
  if size(S,2)<(l+rr+k-1)*n
    S=[S zeros(n,(l+rr+k-1)*n-size(S,2))];
  end; 
  for j=1:k
    for o=-ss:rr
      S(1:n,(l+o-j+k-1)*n+(1:n))=S(1:n,(l+o-j+k-1)*n+(1:n)) + ...
          h*beta(j)*gamma_vect(rr+1-o)*fac*B;
    end;
  end;  
end;

return;
