 function albeck_bayesian_options_model_run(thermo_temp,prior_option,start_file)
cd model_selection

addpath /home/he22/MATLAB/
addpath /home/he22/MATLAB/model_selection
% thermo_temp=1;
% prior_option='general';
% start_file=[];

num_workers = 8;

% create scheduler object
scheduler2 = findResource('scheduler', 'type', 'lsf');
scheduler2.SubmitArguments = '-q sorger_par_unlimited';

% open the pool
matlabpool(scheduler2, num_workers);


global nominal_value_file init_cond_file ode_file prior_flag_file model_ode_observables rel_tol_value likelihood_gen_file

time_max=1200000;
num_nodes=20;
hess_option='hess';
if isempty(start_file)
    start_option='random';
else
    start_option='continue';
end
boundary_option='noboundary';
%prior_option='general';
traj_option=1;
data_file='sabrina_slope2_1_traj.mat';
%thermo_temp=1;
%size_k=[];
simulated_annealing='yes';
first_hessian=25000; %25000;
nominal_value_file=@albeck_nominal_values;
init_cond_file=@albeck_init_conds;
ode_file=@albeck_odes;
prior_flag_file=@albeck_prior_flags;
model_ode_observables = @albeck_observables;
%k=22:116; 
%start_file=[];
rel_tol_value=1e-3;
likelihood_gen_file=[];
model_name='albeck';


bayesian_options_model_selection(time_max, num_nodes, hess_option,start_option, boundary_option,...
    prior_option, traj_option, data_file, thermo_temp, simulated_annealing, first_hessian,...
    nominal_value_file,init_cond_file, ode_file, prior_flag_file, model_ode_observables, model_name,start_file)



 end

