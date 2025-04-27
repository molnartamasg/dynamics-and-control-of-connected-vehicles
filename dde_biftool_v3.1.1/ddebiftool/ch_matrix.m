function Delta=ch_matrix(funcs,xx,par,lambda,varargin)
%% Characteristic matrix and its derivatives
% function Delta=ch_matrx(funcs,x,par,l)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n (either n x (ntau+1), or n x 1)
%	par parameter values
%	lamba complex number at which  charactersitic matrix is computed
%   optional named argument: 'deri', integer (default 0) return derivative
%   of characteristic matrix wrt lambda
% OUTPUT: 
%	D characteristic matrix in C^(n x n)
%
% |xx| is assumed to be equilibrium such that |xx(:,2:end)| are ignored and
% |xx(:,ones(1,(ntau+1))| is used instead.
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
% 
%
%%
default={'deri',0};
options=dde_set_options(default,varargin);
n = size(xx,1); % n = #coordinates,
if ~funcs.tp_del
    taus = par(funcs.sys_tau());
    taus = [0, taus]; % First delay zero
    r=length(taus); % number of delays, r = #delays+1
else % state-dependent delays
    r=funcs.sys_ntau()+1;
    taus = zeros(1,r); % First delay zero
    for i = 2:r
        taus(i) = funcs.sys_tau(i-1,xx(:,ones(1,i)),par);
    end
end
xx=xx(:,ones(1,r));
lfac=[lambda,1,zeros(1,options.deri-1)];
tpow=options.deri;
tfac=(-1)^options.deri;
Delta = lfac(options.deri+1)*eye(n);
for k = 1:r % For every delay
    Ak = funcs.sys_deri(xx,par,k-1,[],[]);
    Delta = Delta -tfac*taus(k)^tpow*Ak*exp(-lambda*taus(k));
end

end
