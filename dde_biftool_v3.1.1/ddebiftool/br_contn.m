function [branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries,varargin)
%% extend DDE-BIFTOOL branch
% function [c_branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries)
% INPUT:
%   funcs problem functions
%	branch initial branch (contains method and initial points)
%	max_tries maximum number of tries
%   optional (named):
%   'plotaxis' (default []): if plotting is on
%   (branch.method.continuation.plot>=1) then the user may specify an axis
%   into which to plot, default [] chooses gca
% OUTPUT:
%	branch extended branch
%	succ number of succesfull corrections
%	fail number of failed corrections
%	rjct number of failures which ended in rejecting a point

% (c) DDE-BIFTOOL v. 3.1.1(125), 29/05/2016
%
% 
%
%% introduce optional argument to set plotting axis 
default={'plotaxis',[]};
options=dde_set_options(default,varargin);
tries=0;
fail=0;
rjct=0;
successful=1;
bound=0;
bound_tau=0;
stop=0;
stop_tau=0;
kontinue=1;
tp_del=funcs.tp_del;

l=length(branch.point);
if l<=1,
    error('br_contn:start','BR_CONTN: could not start branch, length=%d.',l);
end;

method=branch.method.point;
free_par=branch.parameter.free;
max_step=branch.parameter.max_step;
min_bound=branch.parameter.min_bound;
max_bound=branch.parameter.max_bound;

growth_factor=branch.method.continuation.steplength_growth_factor;

% prepare plotting:

if branch.method.continuation.plot>0,
    if branch.method.continuation.plot>=1 && isempty(options.plotaxis)
        options.plotaxis=gca;
    end
end

while kontinue && (tries<=max_tries || bound || bound_tau),
    
    tries=tries+1;
    
    l=length(branch.point);
    
    if l==1
        error('br_contn:fail',['BR_CONTN: could not continue branch,',...
            '%d points, %d fails, % rejected.'],tries-fail,fail,rjct);
    end
    
    prev_point=branch.point(l-1);
    last_point=branch.point(l);
    
    % check boundaries
    
    if bound, % if we are on the boundary, we should stop, unless we crossed more
        stop=1;
        bound=0;
    end;
    if bound_tau, % if delay crossed zero, we should stop
        stop_tau=1;
        bound_tau=0;
    end;
    for j=1:size(min_bound,1),
        if last_point.parameter(min_bound(j,1))<min_bound(j,2), % over minimum
            bound=min_bound(j,1);
            param=last_point.parameter(min_bound(j,1));
            bound_fraction=(param-min_bound(j,2) ) / ...
                (param-prev_point.parameter(min_bound(j,1)));
            bound_parameter=min_bound(j,2);
            break;
        end;
    end;
    for j=1:size(max_bound,1),
        if last_point.parameter(max_bound(j,1))>max_bound(j,2), % over maximum
            bound=max_bound(j,1);
            param=last_point.parameter(max_bound(j,1));
            bound_fraction=(param-max_bound(j,2) ) / ...
                (param-prev_point.parameter(max_bound(j,1)));
            bound_parameter=max_bound(j,2);
            break;
        end;
    end;
    
    if tp_del~=0
        % check sign of delays
        [delay_nr,t_z]=p_tsgn(funcs,last_point);
        if ~isempty(delay_nr)
            bound_tau=1;
        end
    end;
    
    if bound && bound_tau
        bound_tau=0; % we first treat case bound~=0
    end;
    
    if (tries>max_tries) && bound==0 && bound_tau==0
        break;
    end;
    
    if bound, % we already crossed a boundary
        bound_secant=p_axpy(0,last_point,[]);
        bound_secant.parameter(bound)=1;
    elseif (stop && bound_tau==0) || stop_tau
        tries=tries-1;
        break;
    end;
    
    % predict and determine steplength
    
    if l==2 || branch.method.continuation.prediction==1,
        % linear prediction
        secant=p_axpy(-1,last_point,prev_point);
        dist=p_norm(secant);
        if successful, % use extrapolation
            steplength=growth_factor*dist;
        else % use interpolation
            steplength=-dist/2;
        end;
        if bound,
            steplength=-bound_fraction*dist;
        end;
        new_point=p_axpy(-steplength/dist,secant,last_point);
        % check for maximal steplengths
        fraction=1;
        for j=1:size(max_step,1)
            if max_step(j,1)==0
                dp=p_norm(p_axpy(-1,last_point,new_point));
            elseif max_step(j,1)==-1
                if isfield(new_point,'period');
                    lp=setfield(last_point,'period',new_point.period); %#ok<SFLD>
                else
                    lp=last_point;
                end
                dp=p_norm(p_axpy(-1,lp,new_point));
            else
                dp=abs(new_point.parameter(max_step(j,1))-last_point.parameter(max_step(j,1)));
            end
            if dp>max_step(j,2),
                f=max_step(j,2)/dp;
                if f<fraction,
                    fraction=f;
                end;
            end;
        end;
        if fraction<1,
            steplength=steplength*fraction;
            new_point=p_axpy(-steplength/dist,secant,last_point);
        end;
        if bound_tau %negative delay
            % new_point=(tau_p*last_point-tau_n*prev_point)/(tau_p-tau_n)
            tau_n=p_tau(funcs,last_point,delay_nr,t_z);
            tau_p=p_tau(funcs,prev_point,delay_nr,t_z);
            del_tau=tau_p-tau_n;
            pp=p_axpy(-tau_n/del_tau,prev_point,[]);
            new_point=p_axpy(tau_p/del_tau,last_point,pp);
        end;
    else
        err=[branch.method.continuation.prediction];
        error('br_contn:pred',['BR_CONTN: only linear prediction',...
            'is currently implemented, prediction=%d.'],err);
    end;
    %% plot
    online_plot(options.plotaxis,branch,new_point,last_point,j,'g','pred');
    %% correct
    if bound,
        disp('BR_CONTN warning: boundary hit.');
        new_point.parameter(bound)=bound_parameter;
        [new_point,success]=...
            p_correc(funcs,new_point,free_par,bound_secant,method,tries+1,prev_point);
    elseif bound_tau
        s=strcat('BR_CONTN warning: delay number_',num2str(delay_nr),' becomes negative.');
        disp(s);
        [new_point,success]=p_correc(funcs,new_point,free_par,[],method,0,[],delay_nr,t_z);
    else
        if ~branch.method.continuation.steplength_condition,
            secant=[];
        else % normalize secant
            secant=p_secant(secant,p_norm(new_point));
        end;
        [new_point,success]=...
            p_correc(funcs,new_point,free_par,secant,method,tries+1,prev_point);
        
    end
    new_success=success;
    %JS: check for nan, inf
    if isfield(new_point,'profile')
        notisfin=sum(~isfinite(new_point.parameter))+...
            sum(sum(~isfinite(new_point.profile)))+...
            ~isfinite(new_point.period);
        if notisfin>0
            succ=tries-fail;
            return
        end
    end
    if new_success
        % do some normalisations
        new_point=p_normlz(new_point);
        % plot
        online_plot(options.plotaxis,branch,new_point,last_point,j,...
            'b','                        cor');
    else
        fail=fail+1;
    end
    
    % keep new or throw away and maybe throw away last branch point too
    
    if new_success,
        if bound || bound_tau,
            branch.point(l)=new_point;
        elseif successful,
            branch.point(l+1)=new_point;
        else
            branch.point(l+1)=branch.point(l);
            branch.point(l)=new_point;
        end;
    elseif ~successful || bound || bound_tau,
        bound=0;
        bound_tau=0;
        if branch.method.continuation.halt_before_reject==0,
            branch.point=branch.point(1:l-1);
        else
            kontinue=0;
        end;
        rjct=rjct+1;
    end;
    
    successful=new_success;
    
end;

succ=tries-fail;

end
%%
function online_plot(ax,branch,new_point,last_point,step,clr,txt)
mth=branch.method.continuation;
if mth.plot<=0
    return
end
if isempty(mth.plot_measure),
    [x_m,y_m]=df_measr(0,branch);
else
    x_m=mth.plot_measure.x;
    y_m=mth.plot_measure.y;
end
if ~isempty(x_m),
    x1=p_measur(last_point,x_m);
    x2=p_measur(new_point,x_m);
else
    x1=step;
    x2=step+1;
end;
if ~isempty(y_m),
    y1=p_measur(last_point,y_m);
    y2=p_measur(new_point,y_m);
else
    y1=step;
    y2=step+1;
end
if mth.plot>=1
    try
        ish=ishold(ax);
        hold(ax,'on');
        plot(ax,[x1 x2],[y1 y2],clr);
        plot(ax,x2,y2,[clr,'.']);
        if ~ish
            hold(ax,'off');
        end
    catch ME
        warning('br_contn:online_plot',...
            'br_contn: online plotting requested but plot command failed.\n%s',ME.message);
    end
    if mth.plot_progress
        drawnow;
    end
elseif mth.plot_progress
    fprintf('%s: (%g, %g)\n',txt,x2,y2);
end
end