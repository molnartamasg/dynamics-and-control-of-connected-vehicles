function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% (c) DDE-BIFTOOL v. 3.1.1(112), 02/09/2015
%%
if nargin>0
    index=bif_num(biftype,'->');
else
    index=bif_num();
end
end
