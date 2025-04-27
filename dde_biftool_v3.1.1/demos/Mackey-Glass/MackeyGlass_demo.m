%% DDE-Biftool demo Mackey-Glass Equation
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
%
% <html>
% (c) DDE-BIFTOOL v. 3.1.1(73), 31/12/2014
% </html>
%
%% load DDE-Biftool into path
clear
close all
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_psol',...
    '../../ddebiftool_extra_nmfm',...
    '../../ddebiftool_utilities');
%% Enable vectorization
% (disable for speed comparison)
x_vectorize=true;
%% Set user-defined functions
% using |gamma| as constant and (|beta|,|n|,|tau|) as parameters
gamma=1.0;
beta_ind=1;
n_ind=2;
tau_ind=3;
if x_vectorize
    f=@(x,xtau,beta,n)beta*xtau./(1+xtau.^n)-gamma*x;
    funcs=set_funcs(...
        'sys_rhs',@(xx,p)f(xx(1,1,:),xx(1,2,:),p(1),p(2)),...
        'sys_tau',@()tau_ind,...
        'x_vectorized',true);
else
    f=@(x,xtau,beta,n)beta*xtau/(1+xtau^n)-gamma*x; %#ok<UNRCH>
    funcs=set_funcs(...
        'sys_rhs',@(xx,p)f(xx(1,1,:),xx(1,2,:),p(1),p(2)),...
        'sys_tau',@()tau_ind);
end
%% Initial parameters and state
beta0=2;
n0=10;
tau0=0;
x0=(beta0-1)^(1/n0);
%% Initialization of branch of non-trivial equilibria
contpar=tau_ind;
nontriv_eqs=SetupStst(funcs,'x',x0,'parameter',[beta0,n0,tau0],'step',0.1,...
    'contpar',contpar,'max_step',[contpar,0.3],'max_bound',[contpar,10]);
%% Compute and find stability of non-trivial equilibria 
disp('Trivial equilibria');
figure(1);clf
nontriv_eqs=br_contn(funcs,nontriv_eqs,3);
nontriv_eqs=br_stabl(funcs,nontriv_eqs,0,1);
nunst_eqs=GetStability(nontriv_eqs);
ind_hopf=find(nunst_eqs<2,1,'last');
fprintf('Hopf bifurcation near point %d\n',ind_hopf);
%% Continue Hopf bifurcation in two parameters
[hbranch,suc]=SetupHopf(funcs,nontriv_eqs,ind_hopf,...
    'contpar',[beta_ind,tau_ind],'dir',beta_ind,'step',1e-3);
figure(2);clf
hbranch=br_contn(funcs,hbranch,30);
hbranch=br_rvers(hbranch);
hbranch=br_contn(funcs,hbranch,30);
%% Compute L1 coefficient 
% to find if Hopf bifurcation is supercritical (L1<0) or subcritical (L1>0)
[L1,L1low]=HopfLyapunovCoefficients(funcs,hbranch);
fprintf('maximal L1 coefficient along Hopf branch: %g\n',max(L1));
fprintf('max of error estimate for L1 coefficient: %g\n',norm(L1-L1low,'inf'));
%% Branch off at  Hopf bifurcation
disp('Branch off at Hopf bifurcation');
fprintf('Initial correction of periodic orbits at Hopf:\n');
[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',1,'intervals',20,'degree',4,...
    'max_bound',[contpar,20],'max_step',[contpar,0.5]);
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of periodic orbit failed');
end
figure(1);
hold on
per_orb=br_contn(funcs,per_orb,60);
per_orb=br_stabl(funcs,per_orb,0,1);
nunst_per=GetStability(per_orb,'exclude_trivial',true);
%% Find period doubling bifurcations in two parameters
ind_pd=find(diff(nunst_per)==1);
[pdfuncs,pdbranch1,suc]=SetupPeriodDoubling(funcs,per_orb,ind_pd(1),...
    'contpar',[beta_ind,tau_ind],'dir',beta_ind,'step',1e-3);
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of period doubling failed');
end
figure(2);
pdbranch1=br_contn(pdfuncs,pdbranch1,30);
pdbranch=br_rvers(pdbranch1);
pdbranch=br_contn(pdfuncs,pdbranch,30);
%% Check Floquet multipliers 
% (note that Floquet multipliers are often unreliable)
pd1sols=pdfuncs.get_comp(pdbranch.point,'solution');
[nunst_pd,floqpd1,triv_defect,pd1sols]=GetStability(pd1sols,...
    'exclude_trivial',true,'funcs',funcs); %#ok<ASGLU>
fprintf('max defect of Floquet multiplier at -1: %g\n',max(abs(floqpd1+1)));
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
[per2,suc]=DoublePsol(funcs,per_orb,ind_pd(1));
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: branching off at period doubling failed');
end
figure(1);
per2=br_contn(funcs,per2,60);
per2=br_stabl(funcs,per2,0,1);
[nunst_per2,dom,triv_defect]=GetStability(per2,'exclude_trivial',true); 
fprintf('max defect of Floquet multiplier at 1: %g\n',max(triv_defect));
%% Continue period doublings in two parameters for secondary PD
ind_pd2=find(diff(nunst_per2)==1);
[pd2funcs,pdbranch2,suc]=SetupPeriodDoubling(funcs,per2,ind_pd2(1),...
    'contpar',[beta_ind,tau_ind],'dir',beta_ind,'step',1e-3);
if ~suc
    error('MackeyGlassDemo:fail',...
        'MackeyGlassDemo: initialization of 2nd period doubling failed');
end
figure(2);
pdbranch2=br_contn(pdfuncs,pdbranch2,30);
pdbranch2=br_rvers(pdbranch2);
pdbranch2=br_contn(pdfuncs,pdbranch2,30);
%% Check Floquet multipliers along period doubling bifurcation
% (Note that Floquet multipliers are often unreliable.)
pd2sols=pdfuncs.get_comp(pdbranch2.point,'solution');
[nunst_pd2,floqpd2,triv_defect,pd2sols]=GetStability(pd2sols,...
    'exclude_trivial',true,'funcs',funcs);
fprintf('max defect of Floquet multiplier at -1: %g\n',max(abs(floqpd2+1)));

%% Plot of period doubling bifurcation x profiles
bifsols={pd1sols,pd2sols,hbranch.point};
floqpd={floqpd1,floqpd2};
get_par=@(i,k)arrayfun(@(x)x.parameter(i),bifsols{k});
figure(3)
clf;
subplot(3,2,[1,2]);
plot(get_par(beta_ind,1),get_par(tau_ind,1),'.-',...
    get_par(beta_ind,2),get_par(tau_ind,2),'.-',...
    get_par(beta_ind,3),get_par(tau_ind,3),'.-');
legend('PD1','PD2','Hopf');
xlabel('\beta');
ylabel('\tau');
title(sprintf(['Hopf, 1st and 2nd period doubling in Mackey-Glass eqn, ',...
    ' n=%g, gamma=1'],n0));
grid on
for k=1:2
    subplot(3,2,2+k);
    hold on
    for i=1:length(bifsols{k})
        plot(bifsols{k}(i).mesh*bifsols{k}(i).period,bifsols{k}(i).profile,'-');
    end
    hold off
    box on
    grid on
    title(sprintf('PD%d: time profiles of period doubling',k));
    xlabel('t');
    ylabel('x');
    subplot(3,2,4+k);
    plot(1:length(bifsols{k}),floqpd{k}+1,'.-');
    grid on
    title(sprintf('PD%d: dist crit Floq mult from -1',k));
    ylabel('error');
    xlabel('point along branch');
end
%% save
save('MGresults.mat');
