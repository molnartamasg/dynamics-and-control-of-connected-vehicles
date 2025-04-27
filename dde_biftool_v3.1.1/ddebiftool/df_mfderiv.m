%% Approximate higher-order derivative with finite differences
% for normal form calculations in DDE-Biftool
%
% $$ y=\frac{\partial^m}{\partial x^m}
%     f(x_0,\ldots x_{n_\tau},p) X_1 \ldots X_m $$
%
% where |Xk| are |(n x (ntau+1))| matrices.
%
% *WARNING* Numerical approximation of high-order derivatives using finite
% differences is unreliable. Use this function only to test your own
% implementation of |sys_mfderi|.
% 
% (c) DDE-BIFTOOL v. 3.1.1(72), 30/12/2014
%
%%
function [yout,yout2,options]=df_mfderiv(funcs,xx,par,devs,varargin)
%% Input
%
% * funcs: structure containing problem functions (fields sys_rhs and
% x_vectorized are used), can also be only rhs function
% * xx: n x (ntau+1) array in which derivative is taken (1st arg of rhs)
% * par: 1 x np array of system parameters (2nd arg of rhs)
% * devs: cell of length m=|dord|, devs{k} is n x (ntau+1) matrix
% representing kth argument of |dord|-linear form
% * optional arguments, see below
%
%% Output
%
% * yout: n x 1 vector, approximate derivative
% * y2: n x 1 vector, lower-order approximation, use y-y2 as an error
% estimate
% * options: structure mirroring optional inputs used during computations
%% Comments
% The general mixed derivative is computed using the polarization identity,
% reducing it to a sum of higher-order derivatives of single-variable
% functions.
%% Take care of complex deviations recursively
% Avoid calling the user function with complex inputs 
for i=1:length(devs)
    imdevi=imag(devs{i});
    if any(imdevi(:)~=0)
        devs{i}=real(devs{i});
        [y1r,y2r]=df_mfderiv(funcs,xx,par,devs,varargin{:});
        devs{i}=imdevi;
        [y1i,y2i]=df_mfderiv(funcs,xx,par,devs,varargin{:});
        yout=y1r+1i*y1i;
        yout2=y2r+1i*y2i;
        return
    end
end
%% Set parameters
default={...
    'chebyshev_order',24,...       % order of interpolation polynomial (multiple of 4)
    'chebyshev_interval',1e-2,...  % length h of interpolation interval
    'x_vectorized',false,...       % if rhs is vectorized (only used if funcs is rhs)
    'output',1};                   % which output should come first, yout or yout2?                
options=dde_set_options(default,varargin);
dord=length(devs);                            % order of derivative
%% Polarization identity and interpolation at Chebyshev points
% |DF(x) [X_1,X_2,...X_m]= sum(factors(j)*Y_j,j=1..nc)|, where
% mats is a matrix of size m x nc consisting of 0's and 1's, and
% |Y_j=mats(1,j)*X_1+...+mats(m,j)*X_m|.
[mats,factors]=dde_polarization_coeffs(dord);
%% User-defined function
% Check if input is function or |funcs| structure.  If the original user
% function was not vectorized, the wrapper |fvec| (below) is used.
if ~isstruct(funcs) && options.x_vectorized
    rhs=funcs;
elseif ~isstruct(funcs) && ~ options.x_vectorized
    rhs=@(x,p)fvec(funcs,x,p);
elseif ~isfield(funcs,'x_vectorized') || ~funcs.x_vectorized
    rhs=@(x,p)fvec(funcs.sys_rhs,x,p);
else
    rhs=funcs.sys_rhs;
end
%% Create directions for m-order directional derivatives
% Each derivative has the form
%
% $$ \frac{d^m}{d t^m}f(x+t Y_i)$$
%
% where the |Y_i| are called |pdevs| below. They are combinations of the
% user-requested |devs| (determined by |mats|).
ncalls=length(factors);
pdevs=zeros(size(xx,1),size(xx,2),ncalls);
for i=1:ncalls
    for k=1:dord
        pdevs(:,:,i)=pdevs(:,:,i)+mats(k,i)*devs{k};
    end
end
[y,y2]=df_chebnum(rhs,dord,xx,par,pdevs,options.chebyshev_interval,options.chebyshev_order);
fac=factors(:)';
xdim=size(y,1);
ox=ones(1,xdim);
yout=sum(fac(ox,:).*y,2);
yout2=sum(fac(ox,:).*y2,2);
if options.output==2
    tmp=yout;
    yout=yout2;
    yout2=tmp;
end
end
%% Determine coefficients for polarization formula to compute higher-order derivatives
% The combinations in |mats| and |signs| help to reduce an arbitrary-order
% mixed derivative to a linear combination of single-directional
% derivatives.
%%
function [mats,factors]=dde_polarization_coeffs(order)
%% input
%
% * |order|: coefficients for mixed derivatives for degree |n=order| are
% computed
%
%% Output
%
% * |mats|: n x (2^(n-1)-1) matrix P of 1s and 0s (logical)
% * |factors|: 1 x (2^(n-1)-1) array S
%
%% Usage of output
% to obtain multilinear form B(x1,x2,...,xn), sum up S(j)*B(yj,yj,...,yj)
% over j, where yj=P(1)*x1+...P(n)*xn.
%%
mats=logical([1,0]);
factors=[1,-1];
for i=2:order
    nend=size(mats,2);
    mats=[mats,mats;true(1,nend),false(1,nend)]; %#ok<AGROW>
    factors=[factors,-factors]; %#ok<AGROW>
end
mats=mats(:,1:end-1);
factors=factors(1:end-1)/factorial(order);
end
%% pseudo-vectorized version of rhs (used if needed)
function y=fvec(rhs,xx,par)
y1=rhs(xx(:,:,1),par);
nvec=size(xx,3);
y=y1(:,1,ones(nvec,1));
for i=2:nvec
    y(:,1,i)=rhs(xx(:,:,i),par);
end
end