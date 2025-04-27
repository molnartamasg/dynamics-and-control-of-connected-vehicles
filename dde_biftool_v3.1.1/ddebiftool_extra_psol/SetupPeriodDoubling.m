%% SetupPeriodDoubling - Initialize continuation of period doubling bifurcation
%%
function [pdfuncs,pdbranch,suc]=SetupPeriodDoubling(funcs,branch,ind,varargin)
%% 
% Simple wrapper around SetupTorusBifurcation to have a sensible name
% See <SetupTorusBifurcation.html> for description of input and output.
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
% </html>
[pdfuncs,pdbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,varargin{:});
end
