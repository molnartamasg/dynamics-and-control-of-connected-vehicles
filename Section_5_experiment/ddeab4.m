%% 4 step Adams-Bashforth scheme to solve DDEs
function xsolu=ddeab4(fun,delay,history,time)
% dimension of problem
dim=length(history(time(1)));
% calculation of the fixed time step
h=time(2)-time(1);
timehist=fliplr(time(1)-delay:h:time(1));
% initial state containing [x_{0};x_{-1};...;x_{-r}]
z0=history(timehist);
if size(z0,2)~=length(timehist) % repeat initial condition if constant
    z0=repmat(z0,1,length(timehist));
end
z0=z0(:);
% memory allocation for solution
xsolu=zeros(dim,length(time));
% solution of the differential equation
for kk=1:length(time)
    % tnew, xnew: time and solution corresponding to the subsequent timestep
    tnew=time(kk);
    % evaluation of initial conditions
    if kk==1
        znew=z0;
    else
        % evaluation of the right-hand side of the ODE
        xold=zold(1:dim);
        xtauold=zold(end-dim+1:end);
        rhs1=fun(told,xold,xtauold);
        % calculation of the solution
        if kk==2
            xnew = xold+h*rhs1;
        elseif kk==3
            xnew = xold+h*(3*rhs1-rhs2)/2;
        elseif kk==4
            xnew = xold+h*(23*rhs1-16*rhs2+5*rhs3)/12;
        elseif kk>4
            xnew=xold+h*(55*rhs1-59*rhs2+37*rhs3-9*rhs4)/24;
        end
        znew=[xnew;zold(1:end-dim)];
    end
    % store right-hand side expressions for the multi step method
    if kk>3
        rhs4=rhs3;
    end
    if kk>2
        rhs3=rhs2;
    end
    if kk>1
        rhs2=rhs1;
    end
    % store solution
    told=tnew;
    zold=znew;	% state containing [x_{k};x_{k-1};...;x_{k-r}]
    xsolu(:,kk)=znew(1:dim);
end