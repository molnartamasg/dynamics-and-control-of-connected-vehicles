function [J,res]=hopf_jac(funcs,x,omega,v,par,free_par,c)
%% Jacobian for Hopf problem
% function [J,res]=hopf_jac(funcs,x,omega,v,par,free_par,c)
% INPUT:
%   funcs problem function
%	x current Hopf solution guess in R^n
%	omega current Hopf frequency guess in R
%	v current eigenvector guess in C^n
%	par current parameters values in R^p
%	free_par free parameter numbers in N^d 
%	c normalization vector of v in C^(1 x n)
% OUTPUT: 
%	J jacobian in R^(3n+2+s x 3n+1+p)
%	res residual in R^(3n+2+s x 1)

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
n=length(x);
sys_tau=funcs.sys_tau;
sys_rhs=funcs.sys_rhs;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;
sys_dtau=funcs.sys_dtau;

tp_del=funcs.tp_del;
if tp_del==0
  n_tau=sys_tau();
  tau=[0 par(n_tau)];
  m=length(n_tau);
  xx=x(:,ones(m+1,1));
else
  m=sys_ntau();
  xx=x(:,ones(m+1,1));
  t_tau=NaN(1,m);
  for j=1:m
    t_tau(j)=sys_tau(j,xx,par);
  end
  tau=[0 t_tau];
end;

n_fp=length(free_par);

l=1i*omega;

dD=eye(n);
D=l*dD;

for j=0:m
  B=sys_deri(xx,par,j,[],[])*exp(-l*tau(j+1));
  D=D-B;
  dD=dD+tau(j+1)*B;
end;

Dv=D*v;
cv=c*v-1;
dDdxv=zeros(n);
for j=0:m
  for k=0:m
    dDdxv=dDdxv-sys_deri(xx,par,[k j],[],v)*exp(-l*tau(k+1));
    if tp_del~=0 && k>0 
      dDdxv=dDdxv+l*sys_deri(xx,par,k,[],[])*exp(-l*tau(k+1))* ...
            (v*sys_dtau(k,xx,par,j,[]));
    end;
  end;
end;

dDdpv=zeros(n,n_fp);

for k=1:n_fp
  for j=0:m
    dDdp=sys_deri(xx,par,j,free_par(k),[]);
    dDdpv(:,k)=dDdpv(:,k)-dDdp*v*exp(-l*tau(j+1));
    if tp_del==0 && j>0
      if n_tau(j)==free_par(k)
        dDdpv(:,k)=dDdpv(:,k) + ...
		l*sys_deri(xx,par,j,[],[])*v*exp(-l*tau(j+1));
      end;
    elseif tp_del~=0 && j>0 
      d=sys_dtau(j,xx,par,[],free_par(k)); % d=d(tau(j))/d(free_par(k))
          dDdpv(:,k)=dDdpv(:,k) + ...
                l*sys_deri(xx,par,j,[],[])*v*exp(-l*tau(j+1))*d;
    end;
  end;
end;

res(1:n,1)=sys_rhs(xx,par);
res(n+1:2*n,1)=real(Dv);
res(2*n+1:3*n,1)=imag(Dv);
res(3*n+1)=real(cv);
res(3*n+2)=imag(cv);

J=zeros(n,n);
for j=0:m
  J=J+sys_deri(xx,par,j,[],[]);
end;
for j=1:n_fp
  J(1:n,3*n+1+j)=sys_deri(xx,par,[],free_par(j),[]);
end;
J(n+1:2*n,1:n)=real(dDdxv);
J(2*n+1:3*n,1:n)=imag(dDdxv);
J(n+1:2*n,n+1:2*n)=real(D);
J(n+1:2*n,2*n+1:3*n)=-imag(D);
J(2*n+1:3*n,n+1:2*n)=imag(D);
J(2*n+1:3*n,2*n+1:3*n)=real(D);
J(n+1:2*n,3*n+1+(1:n_fp))=real(dDdpv);
J(2*n+1:3*n,3*n+1+(1:n_fp))=imag(dDdpv);
J(3*n+1,n+1:2*n)=real(c);
J(3*n+1,2*n+1:3*n)=-imag(c);
J(3*n+2,n+1:2*n)=imag(c);
J(3*n+2,2*n+1:3*n)=real(c);
J(n+1:2*n,3*n+1)=-real(dD)*imag(v)-imag(dD)*real(v);
J(2*n+1:3*n,3*n+1)=real(dD)*real(v)-imag(dD)*imag(v);

end
