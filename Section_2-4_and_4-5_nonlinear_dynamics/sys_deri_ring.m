%% Jacobians of CCC on a ring
function J=sys_deri_ring(xx,par,nx,np,v,ind,N,CHVs,CAVs,n,Vh,dVh,ddVh,V0,dV0,ddV0,W,dW,dsat)
% free parameters
L=par(ind.L);
alphah=par(ind.alphah);
betah=par(ind.betah);
alpha01=par(ind.alpha01);
beta01=par(ind.beta01);
beta0n=par(ind.beta0n);

% states are stored as
% xx=[v_1(t),v_1(t-sigma),v_1(t-tau);
%     v_2(t),v_2(t-sigma),v_2(t-tau);
%     ...;
%     v_N(t),v_N(t-sigma),v_N(t-tau);
%     h_1(t),h_1(t-sigma),h_1(t-tau);
%     h_2(t),h_2(t-sigma),h_2(t-tau);
%     ...
%     h_{N-1}(t),h_{N-1}(t-sigma),h_{N-1}(t-tau)]

% index of states
hCAVidx = N+CAVs(CAVs<N);
vCAVidx = CAVs;
v1CAVidx = CAVs+1 - (CAVs+1>N)*N;
vnCAVidx = CAVs+n - (CAVs+n>N)*N;
hCHVidx = N+CHVs(CHVs<N);
vCHVidx = CHVs;
v1CHVidx = CHVs+1 - (CHVs+1>N)*N;

% velocities, headways and their derivatives - calculate only when necessary
if ismember(1,nx) || any(ismember([ind.alpha01,ind.beta01,ind.beta0n],np))
    hCAV = xx(hCAVidx,2,:);
    if ismember(N,CAVs)
        hCAV = [hCAV;L - sum(xx([hCAVidx,hCHVidx],2,:),1)];
    end
    vCAV = xx(vCAVidx,2,:);
    v1CAV = xx(v1CAVidx,2,:);
    vnCAV = xx(vnCAVidx,2,:);
elseif ismember(2,nx) || any(ismember([ind.alphah,ind.betah],np))
    hCHV = xx(hCHVidx,3,:);
    if ismember(N,CHVs)
        hCHV = [hCHV;L - sum(xx([hCAVidx,hCHVidx],3,:),1)];
    end
    vCHV = xx(vCHVidx,3,:);
    v1CHV = xx(v1CHVidx,3,:);
end
% hdot = xx([v1CAVidx(CAVs<N),v1CHVidx(CHVs<N)],1,:) - ...
%         xx([vCAVidx(CAVs<N),vCHVidx(CHVs<N)],1,:);

% control inputs - calculate only when necessary
if ismember(1,nx) || any(ismember([ind.alpha01,ind.beta01,ind.beta0n],np)) || (ismember([ind.L],np) && ismember(N,CAVs))
    uCAV=alpha01*(V0(hCAV)-vCAV)+beta01*(W(v1CAV)-vCAV)+beta0n*(W(vnCAV)-vCAV);
elseif ismember(2,nx) || any(ismember([ind.alphah,ind.betah],np)) || (ismember([ind.L],np) && ismember(N,CHVs))
    uCHV=alphah*(Vh(hCHV)-vCHV)+betah*(v1CHV-vCHV);
end

% control input derivatives - calculate only when necessary
% wrt states
if isempty(np) && length(nx)==1
    if ismember(1,nx)
        dudhCAV=alpha01*dV0(hCAV);
        dudvCAV=-(alpha01+beta01+beta0n);
        dudv1CAV=beta01*dW(v1CAV);
        dudvnCAV=beta0n*dW(vnCAV);
    end
    if ismember(2,nx)
        dudhCHV=alphah*dVh(hCHV);
        dudvCHV=-(alphah+betah);
        dudv1CHV=betah;
    end
end
% wrt parameters
if isempty(nx)
    if ismember([ind.L],np)
        if ismember(N,CAVs)
            dudhCAV=alpha01*dV0(hCAV);
        else
            dudhCHV=alphah*dVh(hCHV);
        end
    elseif ismember(ind.alpha01,np)
        dudaCAV=V0(hCAV)-vCAV;
    elseif ismember(ind.beta01,np)
        dudb1CAV=W(v1CAV)-vCAV;
    elseif ismember(ind.beta0n,np)
        dudbnCAV=W(vnCAV)-vCAV;
    elseif ismember(ind.alphah,np)
        dudaCHV=Vh(hCHV)-vCHV;
    elseif ismember(ind.betah,np)
        dudbCHV=v1CHV-vCHV;
    end
end
% wrt states and parameters
if ismember(1,nx)
    if ismember([ind.L],np) && ismember(N,CAVs)
        dudhdhCAV=alpha01*ddV0(hCAV);
    elseif ismember(ind.alpha01,np)
        dudhdaCAV=dV0(hCAV);
        dudvdaCAV=-1;
    elseif ismember(ind.beta01,np)
        dudvdb1CAV=-1;
        dudv1db1CAV=dW(v1CAV);
    elseif ismember(ind.beta0n,np)
        dudvdbnCAV=-1;
        dudvndbnCAV=dW(vnCAV);
    end
elseif ismember(2,nx)
    if ismember([ind.L],np) && ismember(N,CHVs)
        dudhdhCHV=alphah*ddVh(hCHV);
    elseif ismember(ind.alphah,np)
        dudhdaCHV=dVh(hCHV);
        dudvdaCHV=-1;
    elseif ismember(ind.betah,np)
        dudvdbCHV=-1;
        dudv1dbCHV=1;
    end
end
% second derivatives wrt states (using ddW=0)
if length(nx)==2
    if ismember(1,nx)
        dudhdhCAV=alpha01*ddV0(hCAV);
    elseif ismember(2,nx)
        dudhdhCHV=alphah*ddVh(hCHV);
    end
end

% conversion of index pairs to one index
% idxpair=@(idx1,idx2)sub2ind([2*N-1,2*N-1],idx1,idx2);
idxpair=@(idx1,idx2)sub2ind([2*N-1,2*N-1,size(xx,3)],...
    repmat(idx1,1,size(xx,3)),repmat(idx2,1,size(xx,3)),...
    kron(1:size(xx,3),ones(size(idx1))));   % same with vectorization

% Jacobian of system
J=[];
% derivative wrt state
if length(nx)==1 && isempty(np) && isempty(v)
    J=zeros(2*N-1,2*N-1,size(xx,3));
    if nx==0 % derivative wrt x(t)
        J(idxpair(hCAVidx,v1CAVidx(CAVs<N)))=1;
        J(idxpair(hCAVidx,vCAVidx(CAVs<N)))=-1;
        J(idxpair(hCHVidx,v1CHVidx(CHVs<N)))=1;
        J(idxpair(hCHVidx,vCHVidx(CHVs<N)))=-1;
    elseif nx==1 % derivative wrt x(t-sigma)
        if ismember(N,CAVs)
            J(idxpair(vCAVidx(CAVs<N),hCAVidx))=dsat(uCAV(1:end-1,:)).*dudhCAV(1:end-1,:);
            J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCAV(end,:)).*dudhCAV(end,:),1,1,N-1),[1,3,2]);
        else
            J(idxpair(vCAVidx,hCAVidx))=dsat(uCAV).*dudhCAV;
        end
        J(idxpair(vCAVidx,vCAVidx))=dsat(uCAV).*dudvCAV;
        J(idxpair(vCAVidx,v1CAVidx))=dsat(uCAV).*dudv1CAV;
        J(idxpair(vCAVidx,vnCAVidx))=dsat(uCAV).*dudvnCAV;
    elseif nx==2 % derivative wrt x(t-tau)
        if ismember(N,CAVs)
            J(idxpair(vCHVidx,hCHVidx))=dsat(uCHV).*dudhCHV;
        else
            J(idxpair(vCHVidx(CHVs<N),hCHVidx))=dsat(uCHV(1:end-1,:)).*dudhCHV(1:end-1,:);
            J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCHV(end,:)).*dudhCHV(end,:),1,1,N-1),[1,3,2]);
        end
        J(idxpair(vCHVidx,vCHVidx))=dsat(uCHV).*dudvCHV;
        J(idxpair(vCHVidx,v1CHVidx))=dsat(uCHV).*dudv1CHV;
    end
% derivative wrt parameter
elseif isempty(nx) && length(np)==1 && isempty(v)
    J=zeros(2*N-1,1,size(xx,3));
    if np==ind.L            % derivative wrt L
        if ismember(N,CAVs)
            J(N,1,:)=dsat(uCAV(end,:)).*dudhCAV(end,:);
        else
            J(N,1,:)=dsat(uCHV(end,:)).*dudhCHV(end,:);
        end
    elseif np==ind.alphah	% derivative wrt alphah
        J(vCHVidx,1,:)=dsat(uCHV).*dudaCHV;
    elseif np==ind.betah	% derivative wrt betah
        J(vCHVidx,1,:)=dsat(uCHV).*dudbCHV;
    elseif np==ind.alpha01	% derivative wrt alpha01
        J(vCAVidx,1,:)=dsat(uCAV).*dudaCAV;
    elseif np==ind.beta01	% derivative wrt beta01
        J(vCAVidx,1,:)=dsat(uCAV).*dudb1CAV;
    elseif np==ind.beta0n	% derivative wrt beta0n
        J(vCAVidx,1,:)=dsat(uCAV).*dudbnCAV;
%     elseif np==ind.sigma	% derivative wrt sigma
%     elseif np==ind.tau	% derivative wrt tau
    end
% derivative wrt state and parameter (using ddsat=0)
elseif length(nx)==1 && length(np)==1 && isempty(v)
    J=zeros(2*N-1,2*N-1,size(xx,3));
    if nx==1 % derivative wrt x(t-sigma)
        if np==ind.L            % derivative wrt L
            if ismember(N,CAVs)
                J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCAV(end,:)).*dudhdhCAV(end,:),1,1,N-1),[1,3,2]);
            end
        elseif np==ind.alpha01	% derivative wrt alpha01
            if ismember(N,CAVs)
                J(idxpair(vCAVidx(CAVs<N),hCAVidx))=dsat(uCAV(1:end-1,:)).*dudhdaCAV(1:end-1,:);
                J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCAV(end,:)).*dudhdaCAV(end,:),1,1,N-1),[1,3,2]);
            else
                J(idxpair(vCAVidx,hCAVidx))=dsat(uCAV).*dudhdaCAV;
            end
            J(idxpair(vCAVidx,vCAVidx))=dsat(uCAV).*dudvdaCAV;
        elseif np==ind.beta01	% derivative wrt beta01
            J(idxpair(vCAVidx,vCAVidx))=dsat(uCAV).*dudvdb1CAV;
            J(idxpair(vCAVidx,v1CAVidx))=dsat(uCAV).*dudv1db1CAV;
        elseif np==ind.beta0n	% derivative wrt beta0n
            J(idxpair(vCAVidx,vCAVidx))=dsat(uCAV).*dudvdbnCAV;
            J(idxpair(vCAVidx,vnCAVidx))=dsat(uCAV).*dudvndbnCAV;
        end
    elseif nx==2 % derivative wrt x(t-tau)
        if np==ind.L            % derivative wrt L
            if ismember(N,CHVs)
                J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCHV(end,:)).*dudhdhCHV(end,:),1,1,N-1),[1,3,2]);
            end
        elseif np==ind.alphah	% derivative wrt alphah
            if ismember(N,CAVs)
                J(idxpair(vCHVidx,hCHVidx))=dsat(uCHV).*dudhdaCHV;
            else
                J(idxpair(vCHVidx(CHVs<N),hCHVidx))=dsat(uCHV(1:end-1,:)).*dudhdaCHV(1:end-1,:);
                J(N,N+1:2*N-1,:)=-permute(repmat(dsat(uCHV(end,:)).*dudhdaCHV(end,:),1,1,N-1),[1,3,2]);
            end
            J(idxpair(vCHVidx,vCHVidx))=dsat(uCHV).*dudvdaCHV;
        elseif np==ind.betah	% derivative wrt betah
            J(idxpair(vCHVidx,vCHVidx))=dsat(uCHV).*dudvdbCHV;
            J(idxpair(vCHVidx,v1CHVidx))=dsat(uCHV).*dudv1dbCHV;
        end
    end
% second derivative wrt state (using ddsat=0 and ddW=0)
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    J=zeros(2*N-1,2*N-1,size(xx,3));
    if nx(1)==1 % derivative wrt x(t-sigma)
        if nx(2)==1	% derivative wrt x(t-sigma)
            if ismember(N,CAVs)
                J(idxpair(vCAVidx(CAVs<N),hCAVidx))=dsat(uCAV(1:end-1,:)).*dudhdhCAV(1:end-1,:).*squeeze(v(hCAVidx,:,:));
                J(N,N+1:2*N-1,:)=permute(repmat(dsat(uCAV(end,:)).*dudhdhCAV(end,:),1,1,N-1),[1,3,2]).*permute(v(N+1:2*N-1,:,:),[2,1,3]);
            else
                J(idxpair(vCAVidx,hCAVidx))=dsat(uCAV).*dudhdhCAV.*squeeze(v(hCAVidx,:,:));
            end
        end
    elseif nx(1)==2 % derivative wrt x(t-tau)
        if nx(2)==2	% derivative wrt x(t-tau)
            if ismember(N,CAVs)
                J(idxpair(vCHVidx,hCHVidx))=dsat(uCHV).*dudhdhCHV.*squeeze(v(hCHVidx,:,:));
            else
                J(idxpair(vCHVidx(CHVs<N),hCHVidx))=dsat(uCHV(1:end-1,:)).*dudhdhCHV(1:end-1,:).*squeeze(v(hCHVidx,:,:));
                J(N,N+1:2*N-1,:)=permute(repmat(dsat(uCHV(end,:)).*dudhdhCHV(end,:),1,1,N-1),[1,3,2]).*permute(v(N+1:2*N-1,:,:),[2,1,3]);
            end
        end
    end
end
if isempty(J)
    error(['SYS_DERI: requested derivative nx=%d, np=%d, size(v)=%d',...
        'could not be computed!'],nx,np,size(v));
end
end