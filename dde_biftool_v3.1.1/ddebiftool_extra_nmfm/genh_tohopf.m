function hopf=genh_tohopf(funcs,point,freqs) %#ok<INUSL,INUSD>
%% convert generalized Hopf point to Hopf bifurcation point
% function hopf_point=genh_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration (ignored
%   here)
% OUTPUT:
%	hopf_point: hopf point
%
% (c) DDE-BIFTOOL v. 3.1.1(115), 02/09/2015
%
%%
hopf = point;
hopf.kind='hopf';
if isfield(point,'stability');
    hopf.stability=point.stability;
end
end
