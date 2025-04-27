%% Create stability charts for the vehicle chain with ZOH and time delay
% clear; %close all; clc;

%% Parameters
kappa=0.6;
deltat=0.1; % 0.0, 0.1
sigma=0.1;  % 0.0, 0.1, 0.2, 0.3, 0.4, 0.55

% range of parameters
betamin=-2;
betamax=3;
alphamin=-1;
alphamax=4;
betast=200;
alphast=200;
aspect_ratio=1;

% range of frequencies
ommin=0;
ommax=2*pi/deltat/5;
dom=2*pi/deltat/2500;

% grid in plane of parameters
beta_v=linspace(betamin,betamax,betast+1);      % vector of beta values
alpha_v=linspace(alphamin,alphamax,alphast+1);	% vector of alpha values
beta_m = zeros(length(beta_v),length(alpha_v));	% matrix of beta values
alpha_m = zeros(size(beta_m));                  % matrix of alpha values
om_v=ommin:dom:ommax;                           % vector of omega values

% exponential of frequencies
expom=exp(1i*om_v*deltat);

% list of parameters to put on figures
parlist=['\kappa=',num2str(kappa,'%3.2f'),' [1/s]   ',...
         '\Deltat=',num2str(deltat,'%3.2f'),' [s]   ',...
         '\sigma=',num2str(sigma,'%3.2f'),' [s]'];

%% Stability analysis
% resolution of additional delay
r=round(sigma/deltat);

% coefficient matrices
a0=[1,-deltat;0,1];

% construction of parameter-independent parts of system matrices
A=zeros(2*r+4);
A(1:2,1:2)=a0;
A(3:end,1:end-2)=eye(2*r+2);

% matrix of eigenvalues
eigall = zeros(length(beta_v),length(alpha_v),size(A,1));
% matrix of transfer function values
T = zeros(length(beta_v),length(alpha_v),length(om_v));

% computation of eigenvalues and frequency response
for kbeta=1:length(beta_v)
    beta=beta_v(kbeta);
    for kalpha=1:length(alpha_v)
        alpha=alpha_v(kalpha);
        % construction of remaining parts of system matrices
        ar=[-alpha*kappa*deltat^2/2,	 (alpha+beta)*deltat^2/2;
               alpha*kappa*deltat,     -(alpha+beta)*deltat];
        A(1:2,end-1:end)=ar;
        % computation of eigenvalues
        eigall(kbeta,kalpha,:) = eig(A);
        % computation of transfer function values
        T(kbeta,kalpha,:) = ((alpha*kappa*deltat^2/2+beta*deltat)*expom+alpha*kappa*deltat^2/2-beta*deltat)./...
            ((expom.^3-2*expom.^2+expom).*expom.^r...
            +(alpha*kappa*deltat^2/2+(alpha+beta)*deltat)*expom+alpha*kappa*deltat^2/2-(alpha+beta)*deltat);
        % store parameter values
        beta_m(kbeta,kalpha) = beta;
        alpha_m(kbeta,kalpha) = alpha;
    end
    % show progress of computations
    disp(['Computation of row #',num2str(kbeta),'/',...
            num2str(length(beta_v)),' of the stability chart.']);
end

% matrix of critical eigenvalues
eigcr = max(abs(eigall),[],3);
% matrix of maximum transfer function magnitudes
Mmax = max(abs(T(:,:,0<om_v & om_v<2*pi/deltat)),[],3);

%% Plot stability chart
figure(1); clf; hold on; box on;
% plant stability
contourf(beta_m,alpha_m,-eigcr,-[1 1],'r','Linewidth',1.5);
% string stability
contourf(beta_m,alpha_m,-Mmax,-[1 1],'b','Linewidth',1.5);
colormap gray
axis([betamin betamax alphamin alphamax]);
pbaspect([aspect_ratio,1,1]);
xlabel('\beta [1/s]');
ylabel('\alpha [1/s]');
title(['Stability chart of CCC with ZOH and time delay',10,parlist]);
% legend('plant stability','string stability','Location','southeast');

%% Plot eigenvalues and frequency response
% % select parameters of interest
% % figure(1); [beta0,alpha0]=ginput(1);	% by clicking on stability chart
% beta0=0.4; alpha0=0.2;                  % or by giving manually
% 
% % get closest parameter values to the selected ones
% beta0=interp1(beta_v,beta_v,beta0,'nearest');
% alpha0=interp1(alpha_v,alpha_v,alpha0,'nearest');
% 
% % get corresponding eigenvalues and frequency response
% eig_v=squeeze(eigall(beta_v==beta0,alpha_v==alpha0,:));
% M_v=squeeze(abs(T(beta_v==beta0,alpha_v==alpha0,:)));
% 
% % plot eigenvalues
% figure(2); clf; hold on; box on;
% plot(real(eig_v),imag(eig_v),'rx','Linewidth',2);
% plot(cos(0:2*pi/100:2*pi),sin(0:2*pi/100:2*pi),'k--');    % unit circle
% pbaspect([1 1 1]);
% xlabel('Re');
% ylabel('Im');
% title(['Eigenvalue plot of vehicle chain with ZOH and time delay',10,parlist,...
%        '   \beta=',num2str(beta0,'%3.2f'),' [1/s]',...
%        '   \alpha=',num2str(alpha0,'%3.2f'),' [1/s]']);
% 
% % plot frequency response
% figure(3); clf; hold on; box on;
% plot(om_v,M_v,'b');
% plot([ommin ommax],[1 1],'k--');    % plot unit gain
% xlim([ommin ommax]);
% xlabel('\omega [rad/s]');
% ylabel('|T(e^{j\omega\Deltat})|');
% title(['Frequency response of vehicle chain with ZOH and time delay',10,parlist,...
%        '   \beta=',num2str(beta0,'%3.2f'),' [1/s]',...
%        '   \alpha=',num2str(alpha0,'%3.2f'),' [1/s]']);