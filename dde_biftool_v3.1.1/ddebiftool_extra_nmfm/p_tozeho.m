function [zeho, success] = p_tozeho(point)
%%  Convert to zero-Hopf point
% function zeho = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	zeho: uncorrected starting guess for zero hopf point
%   success: whether conversion was successful
%
% (c) DDE-BIFTOOL v. 3.1.1(66), 23/12/2014
%
%%
% Set success
success = 1;

zeho = point;

if strcmp(point.kind, 'hopf')
   zeho.kind = 'zeho';
   zeho.flag = '';
   if ~isfield(zeho,'nmfm')
      zeho.nmfm = [];
   end
else
   fprintf('P_TOZEHO: only hopf points can be converted into zero hopf.\n');
   success = 0;
end
end
