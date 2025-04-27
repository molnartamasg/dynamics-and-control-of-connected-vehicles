function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,points)
%% creates cloud of points indicating area where eigenvalue can be
% function pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0)
% INPUT:
%	AA,tau,real_min,delta_real_min,nb_nu,eigA0
% OUTPUT:
%	pts
%   
% COMMENT: is implemented for arbitrary m>=1 in new version but
% computational cost increases exponentially with m m is number of delays)
%
% for formulas and motivations see 
%
% Verheyden, Luzyanina, Roose: Efficient computation of 
% characteristic roots of delay differential equations using LMS methods,
% Journal of Computational and Applied Mathematics 214 (2008) 209 – 226.
%
%
% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
% Added on 05/03/2007
%
% 
%
  
theta=linspace(0,2*pi,nb_nu+1);
theta=theta(1:end-1);
nu=complex(sin(theta),cos(theta));

r_min_d=real_min-delta_real_min;
ss=exp(-real_min*tau);

%% (JS) call to generalized version map_pts of hlp_get_pts_h_new (both fcns below)
%E=hlp_get_pts_h_new(AA,nu,ss); %original
E=map_pts(AA,nu,ss); 
%%
points=[points,E];

rp=real(points);
ip=imag(points);
inds_p=find(rp>=r_min_d);
rp=rp(inds_p);
ip=ip(inds_p);
%

pts=complex(rp,ip);
%
%% ALTERNATIVE:
%% Note: convhull's option 'Pp' does not always 
%% work in the case that "the initial hull is narrow" ...
%% rp=[rp,rp];
%% ip=[ip,-ip];
%% %
%% K=convhull(rp,ip,{'Qt','Pp'});
%% conv_hull=complex(rp(K),ip(K));
%% %
%% pts=conv_hull(2:end);
%% % Could throw away half of these points actually ...
%

end


function E=map_pts(CC,nu,ss)
%% maps points on unit circle to eigenvalues
m=length(CC)-1;           % number of delays
nb_nu=length(nu);         % number of points on unit circle
count=BaseFromNum(nb_nu,m);
count=count(:,count(end,:)<=ceil((nb_nu+1)/2));
nc=size(count,2);
dim=size(CC{1},1);
E=NaN(dim,nc); % array of approximate Eigenvalue locations
ss=[1,ss(:)'];
Mlist=reshape(cat(3,CC{:}),dim*dim,m+1).*ss(ones(dim*dim,1),:);
nu=nu(:);
for k=1:size(count,2)
    nusel=[1;nu(count(:,k))];
    M=Mlist*nusel;
    E(:,k)=eig(reshape(M,dim,dim));
end
E=E(:).';
end
%%
function count=BaseFromNum(npoints,ndelays)
%% create counts for npoints points and ndelays delays
k=0:npoints^ndelays-1;
count=NaN(ndelays,length(k));
count(1,:)=mod(k,npoints);
for j=2:ndelays
    k=k-count(j-1,:);
    k=k/npoints;
    count(j,:)=mod(k,npoints);
end
count=count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% old version, not used anymore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=hlp_get_pts_h_new(CC,nu,ss) %#ok<DEFNU>
m=length(CC)-1;  
nb_nu=length(nu);

E=[];

% m ...
switch m
 case 1
  for i2=1:ceil((nb_nu+1)/2)
    e=eig(CC{1}+ss(1)*nu(i2)*CC{2});
    E=[E,e.'];
  end
 case 2
  for i2=1:nb_nu
    for i3=1:ceil((nb_nu+1)/2)
      e=eig(CC{1}+ss(1)*nu(i2)*CC{2}+ss(2)*nu(i3)*CC{3});
      E=[E,e.'];
    end
  end        
 case 3
  for i2=1:nb_nu
    for i3=1:nb_nu
      for i4=1:ceil((nb_nu+1)/2)
	e=eig(CC{1}+ss(1)*nu(i2)*CC{2}+ss(2)*nu(i3)*CC{3}+ss(3)*nu(i4)*CC{4});
	E=[E,e.'];
      end
    end
  end     
end

end
