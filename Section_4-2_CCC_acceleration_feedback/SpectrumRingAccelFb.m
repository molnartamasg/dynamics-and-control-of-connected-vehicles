%% Create spectrum plots for the ring configuration with acceleration feedback
clear; %close all; clc;

%% Parameters and system definition
% number of vehicles
NN=24;
% vehicle parameters
beta=0.25;
alpha=0.4;
kappa=0.6;
sigma=0.6;
gamma=0.25; % 0, 0.25, 0.5, 0.75, 1, 1.25

% numerical parameters
deltat=0.01;            % time step
r=round(sigma/deltat);	% delay resolution

% range of plot
Remin=-5;	% -5, -0.8,  -6
Remax=1;	%  1,  0.2,   1
Immin=-3;	% -3,   -2, -15
Immax=3;	%  3,    2,  15

%% Semi-discretization
% system matrices
a0=[0,-1;0,0];
asig=[0,0;alpha*kappa,-(alpha+beta)];
b0=[1;0];
bsig=[0;beta];
bsighat=[0;gamma];
c=[0,1];
IN=eye(NN);
RN=IN; RN=[RN(2:end,:);RN(1,:)];
A0=kron(IN,a0)+kron(RN,b0*c);
Asig=kron(IN,asig)+kron(RN,bsig*c);
Asighat=kron(RN,bsighat*c);
dim=size(A0);

% semidiscretization matrices
P=expm(A0*deltat);
if det(A0)==0
    R00=integral(@(s)-(s-sigma+(r-1)*deltat)/deltat*expm(A0*(deltat-s)),...
        0,deltat,'ArrayValued',true);
    R10=integral(@(s)(s-sigma+r*deltat)/deltat*expm(A0*(deltat-s)),...
        0,deltat,'ArrayValued',true);
else
    invA0=inv(A0);
    R00 = (invA0+1/deltat*(invA0^2-(sigma-(r-1)*deltat)/A0)*(eye(dim)-P));
    R10 = (-invA0+1/deltat*(-invA0^2+(sigma-r*deltat)/A0)*(eye(dim)-P));
end
R0=R00*Asig;
R1=R10*Asig;
R0hat=-R00*Asighat/deltat;
R1hat=(R00-R10)*Asighat/deltat;
R2hat=R10*Asighat/deltat;

% construction of state transition matrix
G=zeros(dim*(r+1));
G(dim+1:end,1:end-dim)=eye(dim*r);
G(1:dim,1:dim)=P;
G(1:dim,r*dim+1:(r+1)*dim)=G(1:dim,r*dim+1:(r+1)*dim)+R0+R0hat;
G(1:dim,(r-1)*dim+1:r*dim)=G(1:dim,(r-1)*dim+1:r*dim)+R1+R1hat;
G(1:dim,(r-2)*dim+1:(r-1)*dim)=G(1:dim,(r-2)*dim+1:(r-1)*dim)+R2hat;
% computation of critical eigenvalue
mu=eig(G);
lambda=log(mu)/deltat;
lambda0=lambda(abs(lambda)<1e-8);                       % zero eigenvalues
lambdastab=lambda(abs(lambda)>=1e-8 & real(lambda)<0);  % stable eigenvalues
lambdaunstab=lambda(abs(lambda)>=1e-8 & real(lambda)>0);% unstable eigenvalues

%% Spectrum plot
figure(1); clf; hold on; box on;
plot([Remin,Remax],[0,0],'k--');
plot([0,0],[Immin,Immax],'k--');
% plot(real(lambda),imag(lambda),'b*');
plot(real(lambda0),imag(lambda0),'b*');
plot(real(lambdastab),imag(lambdastab),'g*');
plot(real(lambdaunstab),imag(lambdaunstab),'r*');
axis([Remin,Remax,Immin,Immax]);
pbaspect([1,1,1]);
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title(['Spectrum of the ring configuration with acceleration feedback',10,...
       '\kappa=',num2str(kappa,'%3.2f'),'   \sigma=',num2str(sigma,'%3.2f'),...
       '   \beta=',num2str(beta,'%3.2f'),'   \alpha=',num2str(alpha,'%3.2f'),...
       '   \gamma=',num2str(gamma,'%3.2f')]);