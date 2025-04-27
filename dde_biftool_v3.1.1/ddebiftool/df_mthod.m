function method=df_mthod(funcs,kind,flag_newhheur)
%% create methed structure with default parameters for continuation, solving and stability
% function method=df_mthod(kind,flag_newhheur)
% INPUT:
%   funcs problem functions
%	kind the kind of default method wanted
%       flag_newhheur (optional, default: 1) boolean: use new steplength heuristic
%                     Only used if kind==stst
% OUTPUT:
%	method default method
% COMMENT:
%       The order of the LMS method used in the computation of the
%       characteristic roots is hard coded in this file.
  
% (c) DDE-BIFTOOL v. 3.1.1(109), 31/08/2015
% Update on 05/03/2007 ("flag_newhheur" <=> (imag(method.stability.lms_parameter_rho)~=0) )
%
% 
%
%%
% during continuation:
method.continuation.steplength_condition=1; % use steplength condition
method.continuation.plot=1; % plot new points
method.continuation.prediction=1; % use secant prediction
method.continuation.steplength_growth_factor=1.2; % grow steplength with factor 1.2
method.continuation.plot_progress=1; % plot progress gradually
method.continuation.plot_measure=[]; % use default plot measures
method.continuation.halt_before_reject=0; % rejection of points allowed

% BW+MB: addition
method.bifurcation.radial_tolerance_factor = 0.25;
method.bifurcation.minimal_real_part = -0.1;
method.bifurcation.correction_tolerance = 1e-6;
method.bifurcation.secant_iterations = 30;
method.bifurcation.secant_tolerance = 1e-6;
method.bifurcation.imagthreshold = 1e-6;
method.bifurcation.monitor_eigenvalues = 0;
method.bifurcation.plot_testfunctions = 0;
method.bifurcation.pause_on_bifurcation= 0;
% End addition

switch kind
    case 'stst'
        if ~exist('flag_newhheur','var'),
            flag_newhheur=1;
        end;
        % stst
        method.point.newton_max_iterations=5;
        method.point.newton_nmon_iterations=1;
        method.point.halting_accuracy=1e-10;
        method.point.minimal_accuracy=1e-8;
        method.point.extra_condition=0;
        method.point.print_residual_info=0;
        % root
        order=4;
        delta_region=0.1;
        if flag_newhheur,
            % maximal order LMS method
            lms_method='mxo';
        else
            % BDF method
            lms_method='bdf';
        end;
        [alpha,beta]=time_lms(lms_method,4);
        method.stability.lms_parameter_alpha=alpha;
        method.stability.lms_parameter_beta=beta;
        if flag_newhheur,
            [a_ell,b_ell]=get_lms_ellipse_new(order,delta_region);
            method.stability.lms_parameter_rho=complex(a_ell,b_ell);
        else
            method.stability.lms_parameter_rho=time_saf(alpha,beta,0.01,0.01);
        end;
        method.stability.interpolation_order=order;
        method.stability.max_number_of_eigenvalues=100;
        method.stability.minimal_time_step=0.01;
        method.stability.maximal_time_step=0.1;
        method.stability.minimal_real_part=[];
        method.stability.max_newton_iterations=6;
        method.stability.root_accuracy=1e-6;
        method.stability.remove_unconverged_roots=1;
        if ~flag_newhheur
            method.stability.newheuristics_tests=0;
        end
        if funcs.tp_del
            method.stability.delay_accuracy=-1e-8;
        end;
    case 'fold'
        % fold
        method.point.newton_max_iterations=5;
        method.point.newton_nmon_iterations=1;
        method.point.halting_accuracy=1e-9;
        method.point.minimal_accuracy=1e-7;
        method.point.extra_condition=0;
        method.point.print_residual_info=0;
        % root
        [alpha,beta]=time_lms('bdf',4);
        method.stability.lms_parameter_alpha=alpha;
        method.stability.lms_parameter_beta=beta;
        method.stability.lms_parameter_rho=time_saf(alpha,beta,0.01,0.01);
        method.stability.interpolation_order=4;
        method.stability.max_number_of_eigenvalues=100;
        method.stability.minimal_time_step=0.01;
        method.stability.maximal_time_step=0.1;
        method.stability.minimal_real_part=[];
        method.stability.max_newton_iterations=6;
        method.stability.root_accuracy=1e-6;
        method.stability.remove_unconverged_roots=1;
        if funcs.tp_del
            method.stability.delay_accuracy=-1e-8;
        end;
    case 'hopf'
        % hopf
        method.point.newton_max_iterations=5;
        method.point.newton_nmon_iterations=1;
        method.point.halting_accuracy=1e-9;
        method.point.minimal_accuracy=1e-7;
        method.point.extra_condition=0;
        method.point.print_residual_info=0;
        % root
        [alpha,beta]=time_lms('bdf',4);
        method.stability.lms_parameter_alpha=alpha;
        method.stability.lms_parameter_beta=beta;
        method.stability.lms_parameter_rho=time_saf(alpha,beta,0.01,0.01);
        method.stability.interpolation_order=4;
        method.stability.max_number_of_eigenvalues=100;
        method.stability.minimal_time_step=0.01;
        method.stability.maximal_time_step=0.1;
        method.stability.minimal_real_part=[];
        method.stability.max_newton_iterations=6;
        method.stability.root_accuracy=1e-6;
        method.stability.remove_unconverged_roots=1;
        if funcs.tp_del
            method.stability.delay_accuracy=-1e-8;
        end;
    case 'psol'
        % point
        method.point.newton_max_iterations=5;
        method.point.newton_nmon_iterations=1;
        method.point.halting_accuracy=1e-8;
        method.point.minimal_accuracy=1e-6;
        method.point.extra_condition=0;
        method.point.print_residual_info=0;
        method.point.phase_condition=1;
        method.point.collocation_parameters=[];
        method.point.adapt_mesh_before_correct=0;
        method.point.adapt_mesh_after_correct=3;
        % mult
        method.stability.max_number_of_eigenvalues=100;
        method.stability.minimal_modulus=0.01;
        method.stability.collocation_parameters=[];
        if funcs.tp_del
            method.stability.delay_accuracy=-1e-8;
        end;
    case 'hcli'
        % point
        method.point.newton_max_iterations=10;
        method.point.newton_nmon_iterations=1;
        method.point.halting_accuracy=1e-8;
        method.point.minimal_accuracy=1e-6;
        method.point.extra_condition=0;
        method.point.print_residual_info=0;
        method.point.phase_condition=1;
        method.point.collocation_parameters=[];
        method.point.adapt_mesh_before_correct=0;
        method.point.adapt_mesh_after_correct=3;
    otherwise
        error('df_mthod" kind %s not recongnized',kind);
end
end

