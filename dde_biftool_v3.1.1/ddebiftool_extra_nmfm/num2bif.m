function biftype = num2bif(index)
%% Convert an index to a bifurcation type
% Return number of supported types on 'count'
%
% (c) DDE-BIFTOOL v. 3.1.1(112), 02/09/2015
%
%%
if nargin>0
    biftype=bif_num(index,'<-');
else
    biftype=bif_num();
end
end
