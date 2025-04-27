function nmfm_printbif_type(str,occurence)
%% shortcut printing out number of occurences of type str
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
if sum(occurence)>0
    fprintf('%d %s ',sum(occurence),str);
end
end
