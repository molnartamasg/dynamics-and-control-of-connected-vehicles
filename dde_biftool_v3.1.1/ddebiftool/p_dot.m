function [y,J]=p_dot(point1,point2,varargin)
%% Compute dot product of two periodic solutions (psol structs)
% 2nd output J is derivative wrt to point1
% by default only the integral between profiles is taken
% point2 is remeshed if meshes are different (may make J invalid)
% optional inputs
%
% derivatives (2d-integer array [n1,n2], default [0,0]): compute scalar
%   product between n1th and n2th derivative of point1 and point2
% ind_comp (array of integers between 1 and size(point1.profile,1), default
%   is 1:size(point1.profile,1)): only include ind_comp components into
%   scalar product
% free_par_ind (array of integers between 1 and length(point1.parameter),
%   default is []): add
%   point1.parameter(free_par_ind)'*point2.parameter(free_par_ind) to scalar
%   product
% period (logical, default is false): include point1.period*point2.period
% into scalar product
%
% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%

%% process options
default={'derivatives',[0,0],'ind_comp',1:size(point1.profile,1),'free_par_ind',[],'period',false};
options=dde_set_options(default,varargin);
if ~strcmp(point1.kind,'psol')|| ~strcmp(point2.kind,'psol')
    error('p_dot: only implemented for psol');
end
%% check degree and mesh
m=point1.degree;
l=(size(point1.profile,2)-1)/m; % number of intervals
if isempty(point1.mesh)
    t=0:1/(l*m):1;
else
    t=point1.mesh;
end
t2=point2.mesh;
%% remesh point2 if necessary
if point2.degree~=m || length(t2)~=length(t) || norm(t2-t)>1e-5
    point2=p_remesh(point2,m,t);
end
%% initialize J
J=point2;
J.profile(:)=0;
J.parameter(:)=0;
J.period=0;
%% compute integral from 0 to 1 for p1^T*p2 and its jacobian
p1=point1.profile;
p2=point2.profile;
y=0;
Jloc=poly_dot(m,options.derivatives);
hpow=1-sum(options.derivatives);
for i=1:l
    t_ind=(i-1)*m+1:i*m+1;
    h=t(t_ind(end))-t(t_ind(1));
    for j=options.ind_comp
        Jfac=h^hpow*p2(j,t_ind)*Jloc';
        rjj=Jfac*p1(j,t_ind)';
        y=y+rjj;
        J.profile(j,t_ind)=J.profile(j,t_ind)+Jfac;
    end
end

%% append dot product of periods
if options.period
    J.period=point2.period;
    y=y+point1.period*point2.period;
end
%% append dot product of parameters if requested
for i=1:length(options.free_par_ind)
    y=y+point1.parameter(options.free_par_ind(i))*point2.parameter(options.free_par_ind(i));
    J.parameter(options.free_par_ind(i))=point2.parameter(options.free_par_ind(i));
end
end

function J=poly_dot(degree,derivatives)
%% calculates matrix J such that int_0^1 p1(t)*p2(t)dt=p1int^T*J*p2int
% where p1int and p2int are given as values of the polynomials on
% equidistant nodes 0,1/degree...1
%
% if derivatives is given (a vector [d1,d2] of integers) then
% the matrix for int_0^1 p1^(d1)(t)*p2^(d2)(t)dt is calculated
if nargin<2
    derivatives=[0,0];
end
J=zeros(degree+1);
iv=inv(vander((0:degree)/degree));%t=(0:degree)/degree;
for i=1:degree+1
    q=zeros(1,degree+1);
    q(i)=1;
    for j=1:degree+1
        p=zeros(1,degree+1);
        p(j)=1;
        c1=p*iv';%polyfit(t,p,degree);
        c2=q*iv';%polyfit(t,q,degree);
        for k=1:derivatives(1)
            c1=polyder(c1);
        end
        for k=1:derivatives(2)
            c2=polyder(c2);
        end
        pprod=conv(c1,c2);
        ppint=polyint(pprod);
        J(j,i)=polyval(ppint,1)-ppint(end);
    end
end
end
