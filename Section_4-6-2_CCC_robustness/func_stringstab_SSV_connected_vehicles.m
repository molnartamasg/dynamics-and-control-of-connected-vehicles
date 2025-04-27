function REIM_D=func_stringstab_SSV_connected_vehicles(ax,par)

REIM_D=zeros(1,size(ax,2));

dtau=par.dtau;
w=par.w;
n=size(par.alpha,1);

HVnum=length(par.humandriver);

%figure
M=zeros(HVnum*4+1,HVnum*4+1,length(w));
m=zeros(5,5,length(w));
bndscomplex=zeros(1,length(w));

for k=1:size(ax,2)
    par.beta(par.betaindex(1),1)=ax(1,k); %beta10
    par.beta(par.betaindex(2),1)=ax(2,k); %beta20
    
    T=func_tfmatrix2(n,w,par.gamma,par.kappa,par.alpha,par.beta,par.xi);
    
    for i=1:length(w)
        detT(i)=det(T(:,:,i));
    end
    
    g=max(abs(detT))-1e-5;
    
    if g>1
        
        REIM_D(1,k)=-1;
        
    else
        for i=1:length(w)
            scale1=diag(repmat([par.dkappa,par.dalpha,par.dbeta,1],1,HVnum));
            Scale=[scale1  zeros(size(scale1,1),1);   
                   zeros(1,size(scale1,2)) 1];

            dVhuman=par.kappah;
            tau=par.tau;
            a=par.alphah; 
            b=par.betah;
            
            s=1i*w(i);
            c=dtau*s./(1+1/3.465*dtau*s);
            %w0=2.363/dtau;
            %c=dtau*s/(1+dtau*s/2)*((s/w0)^2+1.676*(s/w0)+1)/((s/w0)^2+1.37*(s/w0)+1);
            
            D=s^2+(dVhuman*a+s*(a+b))*exp(-s*tau);
            
            m(:,:,i)=1/D*[  -a*exp(-s*tau)          -exp(-s*tau)        -exp(-s*tau)        -1                              s+a*exp(-s*tau);
                s^2+s*b*exp(-s*tau)     -(dVhuman+s)*exp(-s*tau) -(dVhuman+s)*exp(-s*tau) -(s+dVhuman)                         s*dVhuman-s*b*exp(-s*tau);
                -a*s*exp(-s*tau)        -s*exp(-s*tau)      -s*exp(-s*tau)      -s                              s^2+s*a*exp(-s*tau);
                c*s^2*a*exp(-s*tau)     c*exp(-s*tau)*s^2   c*exp(-s*tau)*s^2   c*exp(-s*tau)*(-dVhuman*a-s*(a+b))   c*exp(-s*tau)*(dVhuman*s*a+s^2*b);
                s*a*exp(-s*tau)         s*exp(-s*tau)       s*exp(-s*tau)       s                               exp(-s*tau)*(dVhuman*a+b*s)];
            
            m11=m(1:4,1:4,i);
            m12=m(1:4,5,i);
            m21=m(5,1:4,i);
            if HVnum==1
                M(1:4,1:4,i)=m11;
                M(1:4,5,i)=m12;
                M(5,1:4,i)=m21*T(1,1,i);
                M(5,5,i)=det(T(:,:,i));
            elseif HVnum==2
                M(1:4,1:4,i)=m11;
                M(1:4,5:8,i)=m12*m21;
                M(1:4,9,i)=m12*T(3,3,i);
                M(5:8,5:8,i)=m11;
                M(5:8,9,i)=m12;
                M(9,1:4,i)=m21*T(1,1,i);
                M(9,5:8,i)=m21*(T(2,1,i)+T(1,1,i)*T(2,2,i));
                M(9,9,i)=det(T(:,:,i));
            end
            M(:,:,i)=M(:,:,i)*Scale;
            
            [bndscomplex(i)]=max(eig(abs(M(:,:,i))));
        end
        
        if max(bndscomplex)<1+par.threshold
            musign=1;
        else
            flag=0;
            for i=1:length(w)
                if bndscomplex(i)<1+par.threshold
                elseif flag==1
                else
                    s=1i*w(i);
                    D=s^2+(dVhuman*a+s*(a+b))*exp(-s*tau);
                    if w(i)==0
                        theta=0;
                        %     elseif w(i)*dtau>=pi
                        %         theta=0;
                    else
                        theta=tan(w(i)*dtau/2)/w(i);
                    end
                    scale1=diag(repmat([par.dkappa,par.dalpha,par.dbeta,theta],1,HVnum));
                    Scale=[scale1  zeros(size(scale1,1),1);   
                           zeros(1,size(scale1,2)) 1];

                    m(:,:,i)=1/D*[  -a*exp(-s*tau)          -exp(-s*tau)                -exp(-s*tau)                2                                           s+a*exp(-s*tau);
                        s^2+s*b*exp(-s*tau)     -(dVhuman+s)*exp(-s*tau)    -(dVhuman+s)*exp(-s*tau)    2*(s+dVhuman)                               s*dVhuman-s*b*exp(-s*tau);
                        -a*s*exp(-s*tau)        -s*exp(-s*tau)              -s*exp(-s*tau)              2*s                                         s^2+s*a*exp(-s*tau);
                        s^3*a*exp(-s*tau)       exp(-s*tau)*s^3             exp(-s*tau)*s^3             s*(-s^2+(a*dVhuman+s*(a+b))*exp(-s*tau))    s^2*exp(-s*tau)*(dVhuman*a+s*b);
                        s*a*exp(-s*tau)         s*exp(-s*tau)               s*exp(-s*tau)               -2*s                                        exp(-s*tau)*(dVhuman*a+b*s)];

                    m11=m(1:4,1:4,i);
                    m12=m(1:4,5,i);
                    m21=m(5,1:4,i);
                    if HVnum==1
                        M(1:4,1:4,i)=m11;
                        M(1:4,5,i)=m12;
                        M(5,1:4,i)=m21*T(1,1,i);
                        M(5,5,i)=det(T(:,:,i));
                    elseif HVnum==2
                        M(1:4,1:4,i)=m11;
                        M(1:4,5:8,i)=m12*m21;
                        M(1:4,9,i)=m12*T(3,3,i);
                        M(5:8,5:8,i)=m11;
                        M(5:8,9,i)=m12;
                        M(9,1:4,i)=m21*T(1,1,i);
                        M(9,5:8,i)=m21*(T(2,1,i)+T(1,1,i)*T(2,2,i));
                        M(9,9,i)=det(T(:,:,i));
                    end
                    M(:,:,i)=M(:,:,i)*Scale;
                    
                    deltasetreal=[-ones(size(M,1)-1,1) zeros(size(M,1)-1,1);
                                    1 0];
                    
                    [bnds] = mussv(M(:,:,i),deltasetreal);
                    
                    if min(bnds)>1
                        flag=1;
                        musign=-1;
                    elseif max(bnds)>(1+par.threshold)
                        flag=1;
                        musign=-1;
                    else
                        musign=1;
                    end
                end
            end
        end
        
        REIM_D(1,k)=musign;
    end
end
end

