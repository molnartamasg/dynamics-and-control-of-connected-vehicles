function eqf=auto_eqd(dtm,ups)
%% find cumulative discretization error on given mesh
% function [eqf]=auto_eqd(dtm,ups)
% INPUT:
%   dtm interval lengths of mesh
%   ups solution profile
% OUTPUT:
%       eqf values of monotonically increasing function eqdf
% COMMENT: 
%       this function is a matlab translation of the AUTO
%       fortran function EQDF, except
%	- it is restricted to the case of periodic solutions
%       - it uses different ordering of the data in ups
%	- the integral is approximated using the trapezium rule

% (c) DDE-BIFTOOL v. 3.1.1(50), 11/05/2014
%
% 
%
%% find out dimensions of solution array
%	ntst number of intervals
%   ndim system dimension
%	ncol number of collocation points (per interval)
hmach=1.0d-7;
[ndim,npoints]=size(ups);
ntst=numel(dtm);
ncol=(npoints-1)/ntst;
%% calculate highest derivative of collocation polynomials
% wh: ncol+1 coefficients of central difference formula for ncol derivative
% on unit interval
wh=auto_cnt(ncol);
wh=repmat(wh,ndim,1);

%% hd: highest derivative of collocation polynomials on each subinterval
hd=NaN(ndim,ntst+1); 
for j=1:ntst
    sc=1/(dtm(j)^ncol);
    ibase=(j-1)*ncol;
    hd(:,j)=sum(wh.*ups(:,ibase+(1:ncol+1))*sc,2);
end
%% Take care of "small derivative" case:
if max(abs(hd(:)))<hmach
    eqf=0:ntst;
    return
end
%% extrapolate to get one more interval (exploit periodicity)
hd(:,ntst+1)=hd(:,1);
dtm(ntst+1)=dtm(1);
%% compute approximation to (ncol+1)-st derivative, by taking divided
% differences of neighboring hd values
scav=2./(dtm(1:end-1)+dtm(2:end));
scav=repmat(scav,ndim,1);
hdp1=(hd(:,2:end)-hd(:,1:end-1)).*scav;

%% define the cumulative error err_dt, to be equidistributed
err=sum(abs(hdp1).^(1/(ncol+1)),1);
err_dt=0.5*(err+err([end,1:end-1])).*dtm(1:end-1);
% avoiding zero error
err_dt=max(err_dt,eps);
eqf=[0,cumsum(err_dt)];
end
%% central difference approximation of nth derivative
function d=auto_cnt(n)
% function [d]=auto_cnt(n)
% INPUT:
%	n n-th derivative
% OUTPUT:
%	d coefficients of central difference formula
% COMMENT: 
%       this function is a matlab translation of the AUTO
%       fortran function CNTDIF
%%
d=ones(1,n+1);       
if n==0
    return
end
for i=1:n,
    d(i+1)=nchoosek(n,i);
end
d(end-1:-2:1)=-d(end-1:-2:1);
% Scale to [0,1]  :
d=d*n^n;
end
