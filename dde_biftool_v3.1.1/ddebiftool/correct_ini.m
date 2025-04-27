function [branch,suc]=correct_ini(funcs,branch,pini,dir,step,correc)
%% set and correct first (two) initial points along branch
%
% used in SetupFold, SetupHopf, SetupPOFold, SetupTorusbifurcation,
% SetupPeriodDoubling
%
% [branch,suc]=correct_ini(funcs,branch,pini,dir,step,correc)
%
% inputs
% funcs:  (struct) functions structure defining problems and rhs
% branch: (stuct) branch structure with initially empty point field
% pini:   (struct) initial guess for first point on branch
% dir:    (int) index of parameter value to change for second point 
%         (this routine cannot be used for branching off),
% step:   (double) deviation taken in parameter dir for second point
% correc: (logical) whether to perform Newton correction or not
%
% (c) DDE-BIFTOOL v. 3.1.1(124), 05/02/2016
%
%%
suc=true;
if correc
    %% correct initial guess
    mth=branch.method.point;
    if mth.extra_condition
        rcond=funcs.sys_cond(pini);
        ncond=length(rcond);
    else
        ncond=0;
    end
    if strcmp(pini.kind,'hopf')||strcmp(pini.kind,'fold')
        %% fix for Hopf because add. conditions are not implemented as sys_cond
        ncond=ncond+1;
    end
    [pfirst,suc]=p_correc(funcs,pini,branch.parameter.free(end-ncond+1:end),...
        [],mth,1,pini);
    if ~suc
        pfirst=pini;
        warning('correct_ini:fail','Correction failed');
    end
else
    pfirst=pini;
end
branch=rmfield(branch,'point');
branch.point(1)=pfirst;
if ~suc
    return
end
%% append 2nd point if desired
if ~isempty(dir)
    p2=branch.point(1);
    p2.parameter(dir)=p2.parameter(dir)+step;
    if correc
        corpar=setdiff(branch.parameter.free,dir);
        [p2c,suc]=p_correc(funcs,p2,corpar,[],mth,1,p2);
        if suc
            branch.point(2)=p2c;
        else
            warning('correct_ini:fail','Correction failed');
        end
    else
        p2c=p2;
    end
    branch.point(2)=p2c;
end
end
