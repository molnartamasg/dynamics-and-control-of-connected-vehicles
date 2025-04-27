function REIM_D=func_stringstab(ax,par)

REIM_D=zeros(1,size(ax,2));

for k=1:size(ax,2)
    par.beta(par.betaindex(1),1)=ax(1,k); %beta10
    par.beta(par.betaindex(2),1)=ax(2,k); %beta20

    F=func_tfmatrix2(par.n,par.w,par.gamma,par.kappa,par.alpha,par.beta,par.xi);
    
%      abs(F)
%      pause
    for i=1:size(F,3)
        detF(i)=det(F(:,:,i));
    end
    REIM_D(1,k)=max(abs(detF))-1e-5-1;
       
end

