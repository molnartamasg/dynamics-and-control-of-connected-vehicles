function [res,J]=rot_cond(point,A,usercond)
%% extra phase condition needed to fix rotation phase
%
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
%
userpoint=point;
userpoint.parameter=userpoint.parameter(1:end-1);
[userres,userJ]=usercond(userpoint);
if isfield(point,'x')
    resphas=0;
    Jphas=p_axpy(0,point,[]);
    Jphas.x=A*point.x;
elseif strcmp(point.kind,'psol')
    point.profile=A*point.profile;
    [rdum,J]=p_dot(point,point); %#ok<ASGLU>
    resphas=0;
    Jphas=J;
else
    error('rot_cond:type','rot_cond: type %s not supported',point.kind);
end
res=[userres(:);resphas];
if ~isempty(userJ) % fix for octave
    J=[userJ(:);Jphas];
else
    J=Jphas;
end
end
