function secant=p_secant(secant,norm_point)
%% compute normalized secant
% function n_secant=p_secant(secant,norm_point)
% INPUT: 
%	secant unnormalized secant
%	norm_point norm of a computed reference point
% OUTPUT:
%	n_secant normalized secant

% (c) DDE-BIFTOOL v. 3.1.1(106), 22/08/2015 
%
% 
%
%%
secant=p_axpy(norm_point/p_norm(secant),secant,[]);

if strcmp(secant.kind,'psol') || strcmp(secant.kind,'hcli')
  m=secant.degree;
  ll=(size(secant.profile,2)-1)/m;
  if isempty(secant.mesh)
    mesh=0:1/(ll*m):1;
  else
    mesh=secant.mesh;
  end;
  dmesh(1)=(mesh(2)-mesh(1))/2;
  dmesh(2:ll*m)=(mesh(3:ll*m+1)-mesh(1:ll*m-1))/2;
  dmesh(ll*m+1)=(mesh(ll*m+1)-mesh(ll*m))/2;
  for k=1:size(secant.profile,2)
    secant.profile(:,k)=secant.profile(:,k)*dmesh(k);
  end;
end;

return;
