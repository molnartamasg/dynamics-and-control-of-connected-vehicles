function [psol,stpcond]=p_topsol(funcs,point,ampl,col_degree,nr_int)
%% create starting guess for periodic solution derived from point
% function [psol_point,stpcond]=p_topsol(funcs,point,ampl,col_degree,nr_int)
% INPUT:
%   funcs problem functions
%	point (with stability information if not hcli)
%	ampl amplitude of periodic solution guess
%       % from Hopf:
%	col_degree piecewise polynomial degree for periodic solution
%	nr_int number of intervals for mesh
%       % for period doubling:
%	col_degree collocation parameters (empty for Gauss) 
%       % from connecting orbit:
%       col_degree, nr_int are optional
% OUTPUT:
%	psol_point starting guess for periodic solution derived from point
%	stpcond steplength condition for use in correction of guess

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
psol.kind='psol';
psol.parameter=point.parameter;
stpcond.kind='psol';
stpcond.parameter=0*point.parameter;

switch point.kind,
  case 'hopf',
    mesh=0:1/(col_degree*nr_int):1;
    x=point.x;
    psol.mesh=mesh;
    psol.degree=col_degree;
    stpcond.mesh=mesh;
    stpcond.degree=col_degree;
    v=point.v/norm(point.v);
    for i=1:size(x,1)
      psol.profile(i,:)=x(i)+ampl*(real(v(i))*sin(2*pi*mesh)+imag(v(i))*cos(2*pi*mesh));
      stpcond.profile(i,:)=real(v(i))*sin(2*pi*mesh)+imag(v(i))*cos(2*pi*mesh);
    end;
    if abs(point.omega)>0
      psol.period=abs(2*pi/point.omega);
    else
      disp('P_TOPSOL: zero frequency in Hopf point, period set to zero.');
    end;
    stpcond.period=0;
  case 'psol',
    % remove trivial multiplier
    [i1,i2]=min(abs(point.stability.mu-1)); %#ok<ASGLU>
    point.stability.mu(i2) = 0;
    % find near-bifurcation multiplier
    [i1,i2]=min(abs(abs(point.stability.mu)-1)); %#ok<ASGLU>
    stab_mth=getfield(df_mthod(funcs,'psol'),'stability'); %#ok<GFLD>
    stab_mth.collocation_parameters=col_degree;
    psol=point;
    nremove=1;
    if real(point.stability.mu(i2))<0 && abs(point.stability.mu(i2)+1)<0.3
      [eig_val,eig_vec]=mult_crit(funcs,point,stab_mth,nremove,-1); %#ok<ASGLU>
      eig_vec=[eig_vec,-eig_vec(:,2:end)];
      psol.mesh=[psol.mesh/2,psol.mesh(2:end)/2+0.5];
      psol.profile=[psol.profile,psol.profile(:,2:end)];
      psol.profile=psol.profile+ampl*eig_vec;
      psol.period=2*psol.period;
    elseif real(point.stability.mu(i2))>0 && abs(point.stability.mu(i2)-1)<0.3
      [eig_val,eig_vec]=mult_crit(funcs,point,stab_mth,nremove,+1); %#ok<ASGLU>
      psol.profile=point.profile+ampl*eig_vec;
    else
      error('P_TOPSOL: periodic solution is not close enough to bifurcation.');
    end;
    stpcond.mesh=psol.mesh;
    stpcond.degree=point.degree;
    stpcond.profile=eig_vec;
    stpcond.period=0;
  case 'hcli',
    psol.mesh=point.mesh;
    psol.degree=point.degree;
    psol.profile=point.profile;
    psol.period=point.period;
    stpcond.mesh=point.mesh;
    stpcond.degree=point.degree;
    stpcond.profile=zeros(size(point.profile));
    stpcond.period=1;
  case 'stst', 
    error('P_TOPSOL: stst to psol not supported, convert to hopf first.');
  case 'fold',
    error('P_TOPSOL: fold to psol not supported, convert to hopf first.');
  otherwise,
    err=point.kind;
    error('P_TOPSOL:kind','point kind %s not recognized.',err);
end;

return;

