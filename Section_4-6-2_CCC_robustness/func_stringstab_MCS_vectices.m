function REIM_D=func_stringstab_MCS_vectices(ax,par)

REIM_D=zeros(1,size(ax,2));

w=par.w;
kappa=par.kappa;
alpha=par.alpha;
beta=par.beta;
xi=par.xi;
gamma=par.gamma;
n=par.n;

npar=4;
indices=zeros(npar,2^npar);
Tp=zeros(n,n,length(w),2^npar);

for k=0:npar-1
i=2^k;
indices(k+1,:)=repmat([ones(1,i) -ones(1,i)],1,2^npar/(2*i));
end

for k=1:size(ax,2)
    beta(par.betaindex(1),1)=ax(1,k); %beta10
    beta(par.betaindex(2),1)=ax(2,k); %beta20
    
	Tnom=func_tfmatrix2(n,w,gamma,kappa,alpha,beta,xi);
     
    for i=1:length(w)
        detT(i)=det(Tnom(:,:,i));
    end
    if max(abs(detT))>1+1e-5
        stab=1;
    flag=1;
    else
        stab=-1;
    flag=0;
    end
    
    %for iter=1:nMCS
     for iter=1:size(indices,2)
        if flag==0
            kappap=kappa; alphap=alpha; betap=beta; xip=xi; 
            
            for kk=1:length(par.humandriver)
                ll=par.humandriver(kk)+1;
                kappap(ll,ll)=kappa(ll,ll)+par.dkappa.*indices(1,iter);
                xip(ll,ll)=par.xi(ll,ll)+par.dtau.*indices(2,iter);
                alphap(ll,ll)=par.alpha(ll,ll)+par.dalpha.*indices(3,iter);
                betap(ll,ll)=par.beta(ll,ll)+par.dbeta.*indices(4,iter);
            end

            par.gamma=gamma;
            Tp(:,:,:,iter)=func_tfmatrix2(n,w,gamma,kappap,alphap,betap,xip);
            
            for i=1:length(w)
                detT(i)=det(Tp(:,:,i,iter));
                if abs(detT(i))>1+1e-5
                    stab=1; flag=1;
                    %disp('Instabil')
                %elseif 
                    %stab=-1; flag=0;
                    %disp('Stabil')
                end
            end
        else
        end
    end
    
    REIM_D(1,k)=stab;
    
end

