function t_new=psol_msh(t,m,profile,l_t_new,m_new)
%% create new mesh equidistributing the error
% function [t_new]=psol_msh(t,m,profile,l_t_new,m_new);
% INPUT:
%	t representation points
%	m number of collocation points (per interval)
%	profile solution profile
%	l_t_new size of new mesh (number of intervals)
%	m_new new number of collocation points (per interval)
% OUTPUT:
%	t_new new, adapted mesh

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%% check consistency of mesh t with degree m
l=(length(t)-1)/m;
if l~=floor(l),
    error('PSOL_MSH: t does not contain l=%d intervals of m=%d points!',length(t),m);
end
%% re-distribute coarse mesh
ti=t(1:m:end);
ti_new=auto_msh(profile,ti,l_t_new);
%% insert internal storage points of collocation polynomial
t_new=NaN(1,m_new*l_t_new+1);
for i=1:l_t_new
    t_new(m_new*(i-1)+1:m_new*i)=ti_new(i)+(ti_new(i+1)-ti_new(i))*(0:m_new-1)/m_new;
end
t_new(end)=1;
end
