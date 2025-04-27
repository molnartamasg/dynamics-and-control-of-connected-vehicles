%% Check Floquet multipliers
% (only for orientation, Floquet multipliers are LESS reliable than
% theextended system)
% This script checks the Floquet multipliers along the newly computed
% bifurcating periodic orbits. It uses the field |'get_comp'| of the
% problem-definitino structure for the extended system to extract the
% underlying periodic orbit from the solution of the extended system.
%
% The output is printed, not plotted.
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(20), 11/04/2014
% </html>
%
clear
load('minimal_demo_results');
basebranches={pbranch,trbranch1,trbranch2};
names={'POfold','Torus branch 1','Torus branch 2'};
biffuncs={pfuncs,trfuncs,trfuncs};
for i=1:length(basebranches)
    basebranches{i}.point=biffuncs{i}.get_comp(basebranches{i}.point,'solution');
    basebranches{i}=br_stabl(funcs,basebranches{i},0,1);
end
for i=1:length(basebranches)
    [bifstab{i},crit{i}]=GetStability(basebranches{i},...
        'exclude_trivial',true,'critical',true); %#ok<SAGROW>
    fprintf('max error along %s: %g\n',names{i},max(abs(abs(crit{i})-1)));
end
%%
save('minimal_demo_results.mat','-append','bifstab','crit');
