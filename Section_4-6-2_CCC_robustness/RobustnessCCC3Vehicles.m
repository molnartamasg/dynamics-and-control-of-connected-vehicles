%% Create robust stability charts for CCC with 3 vehicles
clear; clc; % close all;

%% Parameters
% nominal parameters
% HV
kappah=0.6;
alphah=0.1;
betah=0.6;
tau=1.0;
% CAV
kappacav=0.6;
alphacav=0.4;
sigma=0.6;

% axis labels for stability chart
X='\beta_{01} [1/s]';
Y='\beta_{02} [1/s]';

% range of parameters
beta01=linspace(-0.5,1.5,41); 
beta02=linspace(-0.5,1.5,41);
betaindex=[1,2];

% frequency range
w=linspace(0,2*pi,201);        

% connectivity structure
n=2;            % number of follower vehicles
gamma=[1 1      % adjacency matrix (relevant part)
       0 1].';
% parameters according to connectivity structure
kappa=[kappacav 0
       0        kappah].';
alpha=[alphacav	0
       0        alphah].';
beta=[0	0
      0	betah].';
xi=[sigma   sigma
    0       tau].';

% vehicle number of HV
humandriver=[1];%+1 line

% relative uncertainty (w.r.t. nominal parameters, between 0 and 1, *100%)
p=10;   % 0, 5, 10, 15, 20, 25  % adjust par.threshold accordingly
deltakappa=p/100;
deltaalpha=p/100;
deltabeta=p/100;
deltatau=p/100;

% total uncertainty
dkappa=kappah*deltakappa;
dalpha=alphah*deltaalpha; 
dbeta=betah*deltabeta; 
dtau=tau*deltatau; 

% plot parameters and their uncertainty
figure(1); clf;
subplot(8,2,1); hold on; box on;
plot([0 kappah-dkappa kappah-dkappa kappah+dkappa kappah+dkappa 2*kappah],...
     [0 0 1 1 0 0],'k','Linewidth',2);
plot([kappah kappah],[0 1],'-r','Linewidth',2);
xlabel('\kappa_{h}'); yticks([0 1]); axis tight;
subplot(8,2,3); hold on; box on;
plot([0 alphah-dalpha alphah-dalpha alphah+dalpha alphah+dalpha 2*alphah],...
     [0 0 1 1 0 0],'k','Linewidth',2);
plot([alphah alphah],[0 1],'-r','Linewidth',2);
xlabel('\alpha_{h}'); yticks([0 1]); axis tight;
subplot(8,2,5); hold on; box on;
plot([0 betah-dbeta betah-dbeta betah+dbeta betah+dbeta 2*betah],...
     [0 0 1 1 0 0],'k','Linewidth',2);
plot([betah betah],[0 1],'-r','Linewidth',2);
xlabel('\beta_{h}'); yticks([0 1]); axis tight;
subplot(8,2,7); hold on; box on;
plot([0 tau-dtau tau-dtau tau+dtau tau+dtau 2*tau],...
     [0 0 1 1 0 0],'k','Linewidth',2);
plot([tau tau],[0 1],'-r','Linewidth',2);
xlabel('\tau'); yticks([0 1]); axis tight;

% store parameters in one structure
par.kappah=kappah;	par.alphah=alphah;	par.betah=betah;	par.tau=tau;
par.dkappa=dkappa;	par.dalpha=dalpha;	par.dbeta=dbeta;	par.dtau=dtau;
par.kappa=kappa;	par.alpha=alpha;	par.beta=beta;      par.xi=xi;
par.w=w;
par.n=n;
par.gamma=gamma;
par.humandriver=humandriver;
par.betaindex=betaindex;

%% Nominal stability diagram
disp('Nominal stability diagram');
% settings for nominal stability chart calculation
par.ax=[];
par.ax(1).val=beta01;
par.ax(2).val=beta02;
par.Niteration=4;

% calculate nominal stability chart
tic;
bound_fuction_name='func_stringstab';
outputsolnom=func_db_mdbm(par.ax,bound_fuction_name,par.Niteration,par);

% plot nominal stability chart
subplot(8,2,[9 11 13 15]); hold on; box on;
func_db_plot_mdbm(outputsolnom,'k',[1,2]); view([0 0 1]);
xlabel(X); ylabel(Y); legend('Nominal'); drawnow;
toc;

%% Robust stability diagram
disp('Robust stability diagram');
% settings for robust stability chart calculation
par.Niteration=4;
par.threshold=5e-3;   % 5e-3 (p=0, 5, 10, 15) or 1e-6 (p=20, 25)

% calculate robust stability chart
tic;
bound_fuction_name='func_stringstab_SSV_connected_vehicles';
outputsolrob=func_db_mdbm(par.ax,bound_fuction_name,par.Niteration,par);
toc;

% plot robust stability chart
subplot(8,2,[9 11 13 15]); hold on; box on;
func_db_plot_mdbm(outputsolrob,'r',[1,2]); view([0 0 1]);
xlabel(X); ylabel(Y); legend('Nominal','\mu upper bound'); drawnow;

%% Stability chart intersection for different parameters
disp('Intersection stability diagram');
% settings for stability chart intersection calculation
par.Niteration=4;

% calculate stability chart intersection
tic;
bound_fuction_name='func_stringstab_MCS_vectices';
outputsolpmcs=func_db_mdbm(par.ax,bound_fuction_name,par.Niteration,par);

% plot stability chart intersection
subplot(8,2,[9 11 13 15]); hold on; box on;
func_db_plot_mdbm(outputsolpmcs,'b',[1,2]); view([0 0 1]);
xlabel(X); ylabel(Y); legend('Nominal','\mu upper bound','Validation'); drawnow;
toc;

%% Stability charts at the parameter limits only
disp('Perturbed stability diagrams');
% plot other stability boundaries again for reference
subplot(8,2,[2 4 6 8]); hold on; box on;
func_db_plot_mdbm(outputsolnom,'k',[1,2]);view([0 0 1]); drawnow;
try 
    func_db_plot_mdbm(outputsolrob,'r',[1,2]);view([0 0 1]); drawnow;
    func_db_plot_mdbm(outputsolpmcs,'b',[1,2]);view([0 0 1]); drawnow;
catch
end

% settings for stability chart calculation of perturbed systems
par.Niteration=4;

% calculate and plot stability charts for perturbed systems
bound_fuction_name='func_stringstab';
gr=[-1,1];
for i1=1:length(gr)
for i2=1:length(gr)
for i3=1:length(gr)
for i4=1:length(gr)
    % perturbed parameters
    for kk=1:length(par.humandriver)
        ll=par.humandriver(kk)+1;
        par.kappa(ll,ll)=par.kappa(ll,ll)+dkappa.*gr(i1);
        par.xi(ll,ll)=par.xi(ll,ll)+dtau.*gr(i2);
        par.alpha(ll,ll)=par.alpha(ll,ll)+dalpha.*gr(i3);
        par.beta(ll,ll)=par.beta(ll,ll)+dbeta.*gr(i4);
    end

    % stability chart
    tic;
    outputsolp=func_db_mdbm(par.ax,bound_fuction_name,par.Niteration,par);
    func_db_plot_mdbm(outputsolp,'g',[1,2]);view([0 0 1]); drawnow;
    xlabel(X); ylabel(Y); drawnow;
    toc;
    
    % original nominal parameters
    par.kappa=kappa;
    par.xi=xi;
    par.alpha=alpha;
    par.beta=beta;
end
end
end
end

%% Plot mu curve for a selected point
% pick point
subplot(8,2,[2 4 6 8]);
% [beta10t,beta30t]=ginput(1);
beta10t=0.5; beta30t=0.5;
par.beta(par.betaindex(1),1)=beta10t;
par.beta(par.betaindex(2),1)=beta30t;

% show point
try delete(p1)
    delete(t1)
    p1=plot(beta10t,beta30t,'or','MarkerFaceColor','r');
    t1=text(beta10t,beta30t,'P'); drawnow;
catch
    p1=plot(beta10t,beta30t,'or','MarkerFaceColor','r');
    t1=text(beta10t,beta30t,'P'); drawnow;
end
legend('Nominal','\mu upper bound','Validation','P');

% get corresponding transfer function
w2=w;
T=func_tfmatrix2(n,w2,par.gamma,par.kappa,par.alpha,par.beta,par.xi);

HVnum=length(par.humandriver);
% get corresponding mu curve
M=zeros(HVnum*4+1,HVnum*4+1,length(w2));
m=zeros(5,5,length(w2));
g=zeros(1,length(w2));
for i=1:length(w2)
    dVhuman=par.kappah;
    tau=par.tau;
    a=par.alphah;
    b=par.betah;
    
    s=1i*w2(i);
    D=s^2+(dVhuman*a+s*(a+b))*exp(-s*tau);
    
    if w2(i)==0
        theta=0;
        %     elseif w2(i)*dtau>=pi
        %         theta=0;
    else
        theta=tan(w2(i)*dtau/2)/w2(i);
    end
    
	scale1=diag(repmat([par.dkappa,par.dalpha,par.dbeta,theta],1,HVnum));
    Scale=[scale1                   zeros(size(scale1,1),1);   
           zeros(1,size(scale1,2))	1];

	m(:,:,i)=1/D*[-a*exp(-s*tau)          -exp(-s*tau)              -exp(-s*tau)                2                                           s+a*exp(-s*tau);
                  s^2+s*b*exp(-s*tau)     -(kappah+s)*exp(-s*tau)	-(kappah+s)*exp(-s*tau)     2*(s+kappah)                                s*kappah-s*b*exp(-s*tau);
                  -a*s*exp(-s*tau)        -s*exp(-s*tau)            -s*exp(-s*tau)              2*s                                         s^2+s*a*exp(-s*tau);
                  s^3*a*exp(-s*tau)       exp(-s*tau)*s^3           exp(-s*tau)*s^3             s*(-s^2+(a*kappah+s*(a+b))*exp(-s*tau))     s^2*exp(-s*tau)*(kappah*a+s*b);
                  s*a*exp(-s*tau)         s*exp(-s*tau)             s*exp(-s*tau)               -2*s                                        exp(-s*tau)*(kappah*a+b*s)];

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
    
    g(i)=det(T(:,:,i));
end

deltasetreal=[-ones(size(M,1)-1,1) zeros(size(M,1)-1,1);
              1 0];

tic;
[bnds] = mussv(M,deltasetreal); %!!!!!!!!!!!!!!!!!!!!!!!!
toc;

% plot frequency response and mu curve
subplot(8,2,[10 12 14 16]); hold off;
plot(w2,squeeze(bnds(1,1,:)),'-r'); hold on;
plot(w2,squeeze(bnds(1,2,:)),'-m');
plot(w2,abs(g),'-k');
xlabel('\omega [rad/s]'); ylabel('\mu estimations at point P')
plot([0 w2(end)],[ 1 1],'--k');
xlim([0 w2(end)]);
legend('\mu upper bound','\mu lower bound','abs(G_{n,0})','1'); drawnow;

%% Finalize plot
% title to put on figure
problem='CCC with 3 vehicles';
% list of parameters to put on figure
parlist=['HV: \kappa_{\rm h}=',num2str(kappah,'%3.2f'),' [1/s]   ',...
         '\alpha_{\rm h}=',num2str(alphah,'%3.2f'),' [1/s]   ',...
         '\beta_{\rm h}=',num2str(betah,'%3.2f'),' [1/s]   ',...
         '\tau=',num2str(tau,'%3.2f'),' [s]   ',10,...
         'CAV: \kappa_{0}=',num2str(kappacav,'%3.2f'),' [1/s]   ',...
         '\alpha_{01}=',num2str(alphacav,'%3.2f'),' [1/s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];
suptitle(['Robust stability charts of ',problem,' with ',num2str(p),'% uncertainty']);