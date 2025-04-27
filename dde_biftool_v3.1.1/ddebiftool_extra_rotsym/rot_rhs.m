function y=rot_rhs(xx,p,A,expA,user_rhs,user_tau,isvec)
%% right-hand side in rotating coordinates
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
omega=p(end);
userpar=p(1:end-1);
xxrot=xx;
tau_ind=user_tau();
if ~isvec
    for i=2:size(xx,2)
        xxrot(:,i)=expA(-omega*p(tau_ind(i-1)))*xxrot(:,i);
    end
else
    dim=size(xxrot,1);
    nvec=size(xxrot,3);
    for i=2:size(xx,2)
        xxrot(:,i,:)=expA(-omega*p(tau_ind(i-1)))*reshape(xxrot(:,i,:),dim,nvec);
    end
end
y0=user_rhs(xxrot,userpar);
y=y0-reshape(A*omega*reshape(xxrot(:,1,:),[dim,nvec]),[dim,1,nvec]);
end


