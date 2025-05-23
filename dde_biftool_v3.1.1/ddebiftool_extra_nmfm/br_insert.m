%% Detect bifurcation points, compute their normal form and insert into branch
% (with refinement)
%%
function [bifpoints,biflow,branch,indices]=br_insert(funcs,branch,detect,varargin)
%% Input  
% * funcs: probrlem functions
% * branch: branch along which bifurcation points should be identified
% * detect: cell array of same length as branch.point: if entry ind is
% non-empty then a bifurcation between point(ind) and point(ind+1) had been
% provisionally detected. The entry of the cell array is a function of the
% format
% |[bifpoint,bifp2,branch,bifindex]=detect{ind}(funcs,branch,ind+(0:1),...)|.
% It should return the bifurcation point in bifpoint (in the format used
% for these points), a copy bifp2 (possibly with lower-order numerical finite
% differences for normal forms) to estimate errors, the same |branch| now
% refined, as obtained often by bisection or Newton iterations for
% the special point.  The entries in branch.point should all have the
% format as the original branch. |bifindex| is the number in the refined
% |branch| at which  the bifurcation occurs (should be in the same place as
% |bifpoint|).
%
% Examples of prepared functions for detect cell array (in ddebiftool_utilities):
% TakensBogdanovNormalform, ZeroHopfNormalform, CuspNormalform,
% HopfHopfNormalform, GeneralizedHopfNormalform.
%
%% Output
%
% * bifpoints: cell array of bifurcation points (these can have different
% formats)
% * biflow: same as bifpoints, possibly lower order of finite difference
% (if used)
% * branch: refined branch, including (copies of) bifpoints in the format
% of the points in the input branch.point array
% * indices: integer array of same length as bifpoints: entry bifpoints{k}
% equals branch.point(indices(k)).
%
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
%
%% 
default={'print',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
isbif=@(d)cellfun(@(x)~isempty(x),d);
nextbif=@(d,ind)find(isbif(d(ind+1:end)),1,'first')+ind;
nbifs=sum(isbif(detect));
bifpoints={};
biflow={};
indices=[];
curbif=0;
curind=0;
npoints=length(branch.point);
while curbif<nbifs
    curind=nextbif(detect,curind);
    curbif=curbif+1;
    [bifpoints{curbif},biflow{curbif},branch,indices(curbif)]=detect{curind}(...
        funcs,branch,curind+(0:1),pass_on{:},'print',options.print-1); %#ok<AGROW>
    if options.print>0 && ~isempty(bifpoints{curbif}) && ...
            isstruct(bifpoints{curbif}) && isfield(bifpoints{curbif},'kind');
        fprintf('br_insert: detected %d of %d: %s. Normalform:\n',curbif,nbifs,bifpoints{curbif}.kind);
        disp(bifpoints{curbif}.nmfm);
    end
    nnewpoints=length(branch.point)-npoints;
    detect=[detect(1:curind),cell(1,nnewpoints),detect(curind+1:end)];
    curind=curind+nnewpoints;
    npoints=length(branch.point);
end
end
