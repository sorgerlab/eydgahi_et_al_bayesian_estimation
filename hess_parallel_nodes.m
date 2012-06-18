function total_hess=hess_parallel_nodes(position,num_nodes,traj_option,prior_option,thermo_temp,k,data_file)
s=size(k);
s=s(2);

global model_ode_observables ode_file prior_flag_file init_cond_file likelihood_gen_file


try
%    if num_nodes==0
%         jm=findResource('scheduler','type','local');
%         num_nodes=1;
%         job=createJob(jm);
%     else
    if num_nodes>1
        jm = findResource('scheduler','type','lsf');  %Request an LSF job scheduler
        %jm = schedulerOrchestra; % request a DCT-workalike scheduler for LSF (jeremy's code)
        set(jm,'ClusterMatlabRoot','/opt/matlab')  ;   %Specify /opt/matlab as the MATLAB Root dir
        set(jm, 'SubmitArguments', '-q sysbio_15m');
        job=createJob(jm);
        %job.auto_retry = 4;
    end

catch
   disp('Parallel scheduler could not be started. Turn off the')
   disp('paramater RunOnCluster if your software does not support paralell computing.')
   return
end


total_elements=sum(1:s);
row_number=zeros(1,s);
count=s+1-(1:s);
row_number(1)=1;
for i=2:s
    row_number(i)=row_number(i-1)+count(i-1);
end

if num_nodes <= 1
	num_nodes = 1;
	results = cell(1,13);
	% calculate hessian right here
	[results{:}] = hessian_parallel_short(1,total_elements,1,position,row_number,traj_option,prior_option,thermo_temp,k, model_ode_observables, ode_file, prior_flag_file, init_cond_file, data_file, likelihood_gen_file);
elseif mod(total_elements,num_nodes)==0
    num_elements_in_jobs=total_elements/num_nodes;
    start=1;
    finish=num_elements_in_jobs;
    for job_number=1:num_nodes
        output_hess(job_number)=createTask(job,@hessian_parallel_short,13,{start,finish,job_number,position,row_number,traj_option,prior_option,thermo_temp,k, model_ode_observables, ode_file, prior_flag_file, init_cond_file, data_file, likelihood_gen_file});
        start=start+num_elements_in_jobs;
        finish=finish+num_elements_in_jobs;
    end
else
    num_elements_in_last_job=mod(total_elements,num_nodes-1);
    num_elements_in_jobs=(total_elements-num_elements_in_last_job)/(num_nodes-1);
    
    start=1;
    finish=num_elements_in_jobs;
    for job_number=1:num_nodes-1
        %[start,finish]
        output_hess(job_number)=createTask(job,@hessian_parallel_short,13,{start,finish,job_number,position,row_number,traj_option,prior_option,thermo_temp,k, model_ode_observables, ode_file, prior_flag_file, init_cond_file,data_file, likelihood_gen_file});
        start=start+num_elements_in_jobs;
        finish=finish+num_elements_in_jobs;
        %start_total(job_number)=start;
        % finish_total(job_number)=finish;
    end
    job_number=num_nodes;
    start=total_elements-num_elements_in_last_job+1;
    finish=total_elements;
    output_hess(job_number)=createTask(job,@hessian_parallel_short,13,{start,finish,job_number,position,row_number,traj_option,prior_option,thermo_temp,k,model_ode_observables, ode_file, prior_flag_file, init_cond_file,data_file, likelihood_gen_file});
end

if num_nodes > 1
    submit(job);
    while ~waitForState(job, 'finished', 60)
       disp(datestr(now, 31));
       for i=1:job_number
           fprintf(1, '%3d', i);
       end
       fprintf(1, '\n');
       for i=1:job_number
           fprintf(1, '%3s', upper(output_hess(i).State(1)));
       end
       fprintf(1, '\n\n');
    end
    results=getAllOutputArguments(job);
    destroy(job)
end

save('zzz_divided_short.mat')
empty_flag=0;
for i=1:num_nodes
    if size(results{i,1},2)==0
        empty_flag=1;
    end
    results_matrix(i,:)=horzcat(results{i,1},zeros(1,total_elements-size(results{i,1},2)));
    i
end
elements_in_hessian(1:total_elements)=sum(results_matrix(:,1:total_elements),1);

oneD_hessian2=vertcat(results{:,3});
if empty_flag==0
    total_hess=zeros(s,s);
    element_index=1;
    for i=1:s
        for j=i:s
            total_hess(i,j)=elements_in_hessian(element_index);
            total_hess(j,i)=total_hess(i,j);
            element_index=element_index+1;
        end
    end
else
    total_hess=[];
end
save('zzz_divided_short.mat')


end

