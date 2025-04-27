function n=p_norm(point)
%% norm of point for distance measuring
% function n=p_norm(point)
% INPUT:
%	point 
% OUTPUT:
%	n norm of point

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
% 
%%
switch point.kind,
  case 'stst',
    n=sqrt(norm(point.parameter)^2+norm(point.x)^2);
  case 'fold',
    n=sqrt(norm(point.parameter)^2+norm(point.x)^2+norm(point.v)^2);
  case 'hopf',
    n=sqrt(norm(point.parameter)^2+norm(point.x)^2+norm(point.v)^2+...
      norm(point.omega)^2);  
  case 'psol',
    n=sqrt(norm(point.parameter)^2+...
      norm(point.profile)^2/size(point.profile,2) ...
      +norm(point.period)^2);
  case 'hcli',
    n=sqrt(norm(point.parameter)^2+norm(point.profile)^2/...
      size(point.profile,2)...
      +norm(point.period)^2+norm(point.x1)^2+norm(point.x2)^2+...
      norm(point.lambda_v)^2+norm(point.lambda_w)^2+...
      norm(point.v)^2/size(point.v,2)+norm(point.alpha)^2);
  otherwise,
    err=point.kind;
    error('P_NORM: point %s is not recognized.',err);
end;

end
