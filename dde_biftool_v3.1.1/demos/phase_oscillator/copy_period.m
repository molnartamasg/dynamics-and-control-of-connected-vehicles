%% Use period as parameter (which can be fixed or released)
%
% (c) DDE-BIFTOOL v. 3.1.1(126), 05/09/2016
%
%%
function [res,J]=copy_period(pt,indpar)
res=pt.period-pt.parameter(indpar);
J=p_axpy(0,pt,[]);
J.period=1;
J.parameter(indpar)=-1;
end