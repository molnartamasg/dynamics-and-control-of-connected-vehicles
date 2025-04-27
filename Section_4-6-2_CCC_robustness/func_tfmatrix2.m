function T=func_tfmatrix2(n,w,gamma,kappa,alpha,beta,xi)
T=zeros(n,n,length(w));
s=1i*w;

for i=1:n-1
    T(i,i+1,:)=-1;
end

for j=1:n
    for i=0:j-1
        phi(j,i+1)=alpha(j,i+1)*kappa(j,i+1)/(j-i);
        kappaadd(j,i+1)=alpha(j,i+1)+beta(j,i+1);
    end
end

for j=1:n
    for i=1:j
        
        dsum=0;
    for k=i:n
        dsum=dsum+gamma(k,i)*(s*(alpha(k,i)+beta(k,i))+phi(k,i)).*exp(-s*xi(k,i)); 
    end
        T(j,i,:)=gamma(j,i)*(s*beta(j,i)+phi(j,i)).*exp(-s*xi(j,i))./(s.^2+dsum);
    end
end
