function bayesian_options_model_selection(time_max, num_nodes, hess_option,start_option, boundary_option,...
    prior_option, traj_option, data_file, thermo_temp, simulated_annealing, first_hessian,...
    nominal_value_file,init_cond_file, ode_file, prior_flag_file, model_ode_observables,...
    model_name, start_file)

%Bayesian estimation algorithm used to optimize S. Gaudet's TRAIL model
%--------------------------------------------------------------------------
tic
%global init_cond_file ode_file nominal_value_file prior_flag_file model_ode_observables

if isempty(thermo_temp)
    thermo_temp=1;
end

param=nominal_value_file();
num_param=size(param,2);
num_ic=size(find(prior_flag_file()==0),2);
k=find(prior_flag_file()~=0);
size_k=size(k,2);

switch lower(start_option)
    case 'continue'
        load(start_file,'num_nodes','hess_option','start_option','boundary_option',...
            'prior_option','traj_option','data_file','thermo_temp','size_k');
        start_option='continue';
end

if getenv('LSB_JOBID')
  filename = strcat('log-lsf-', datestr(now, 30), '-', getenv('LSB_JOBID'),'.txt');
else
  filename = strcat('log-', datestr(now, 30) ,'.txt');
end
fileid = fopen(filename, 'w');
fprintf(fileid, '%s p=%d\n', time_max,num_nodes,hess_option,start_option,boundary_option,prior_option,traj_option);

start_time=floor(sum(clock*1000));
rand('state',start_time);
randn('state',start_time);

norm_step_size=0.75;

acceptance=0;
% nominal_vals=nominal_value_file(); %[kc, kf, kr, v, kdeg] = values;

% if traj_option == 1
%         [data_t, data, mito_start, mito_end]=data_file();
% elseif or(traj_option==2, traj_option==3)
%         [data_t, data_ec,mito_start,mito_end]=sabrina_slope_EC;
%         [data_t, data_ic]=sabrina_slope_IC;
%         data=[data_ec,data_ic];
% end

%data_t=data_t';
sig_value=1;%2.63;
T=50;  %the lower the T, the lower the acceptance rate.
T_init=50;

%-------------INITIALIZIATION-------------%
old_posterior=zeros(1,time_max);
rej_post=zeros(1,time_max);
d_post=zeros(1,time_max);
priors=zeros(1,time_max);
lkls=zeros(1,time_max);
alpha_matrix=zeros(1,time_max);
test=zeros(1,time_max);
sigma_matrix=zeros(1,time_max);
delta_lkls=zeros(1,time_max);
delta_post_matrix=zeros(1,time_max);
delta_prior=zeros(1,time_max);
lkls_pre_thermo=zeros(1,time_max);
new_pos=zeros(time_max,num_param);
reject=zeros(1,time_max);
test_lkls=zeros(1,time_max);
test_position=zeros(time_max,num_param);
test_priors=zeros(1,time_max);
total_time=zeros(1,time_max);
%-----------------------------------------%

switch lower(start_option)
    case 'random'
        mult_factor=zeros(1,num_param);
        mult_factor(k)=unifrnd(-1,1,[1,size_k]);
        mult_factor=10.^mult_factor;
        position=mult_factor.*nominal_value_file();%position=actual_init_data(mult_factor);
        position=log10(position);
        disp('random start')
        start_p=1;
        if and(thermo_temp<1,thermo_temp>0)
            if traj_option>1
                thermo_temp_str=num2str(thermo_temp);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_0',thermo_temp_str(3),'_',getenv('LSB_JOBID'),'.mat'];
            elseif traj_option==1
                file_name=char(data_file);
                thermo_temp_str=num2str(thermo_temp);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',file_name(1),file_name(end),'_0',thermo_temp_str(3),'_',getenv('LSB_JOBID'),'.mat'];
            end
        elseif or(thermo_temp==1,thermo_temp==0)
            if traj_option>1
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',num2str(thermo_temp),'_',getenv('LSB_JOBID'),'.mat'];
            elseif traj_option==1
                file_name=char(data_file);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',file_name(1),file_name(end),'_',num2str(thermo_temp(end)),'_',getenv('LSB_JOBID'),'.mat'];
            end
        end
        
        
        
    case 'same'
        %load('10krun_samestart_boundfixed_1.mat', 'new_pos');
        %position=new_pos(9000,:);
        %clear new_pos;
        %load('hessian_init_calc.mat','final_hess_c');
        
        position=nominal_value_file();% position=actual_init_data(ones(1,size_k));
        position=log10(position);
        disp('same start')
        start_p=1;
        if and(thermo_temp<1,thermo_temp>0)
            if traj_option>1
                thermo_temp_str=num2str(thermo_temp);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_0',thermo_temp_str(3),'_',getenv('LSB_JOBID'),'.mat'];
            elseif traj_option==1
                file_name=char(data_file);
                thermo_temp_str=num2str(thermo_temp);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',file_name(1),file_name(end),'_0',thermo_temp_str(3),'_',getenv('LSB_JOBID'),'.mat'];
            end
        elseif or(thermo_temp==1,thermo_temp==0)
            if traj_option>1
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',num2str(thermo_temp),'_',getenv('LSB_JOBID'),'.mat'];
            elseif traj_option==1
                file_name=char(data_file);
                output=[model_name,'_', num2str(traj_option),'_',hess_option(1),start_option(1),boundary_option(1),'_',prior_option(1:3),'_',file_name(1),file_name(end),'_',num2str(thermo_temp(end)),'_',getenv('LSB_JOBID'),'.mat'];
            end
        end
    case 'continue'
        adjusted_time_max=time_max;
        load(start_file);
        start_p=size(new_pos,1);
        position=new_pos(start_p,:);
        disp('continue start')
        clear time_max;
        time_max=adjusted_time_max;
end

init_pos=position;

old_prior=prior_gen(position, prior_option,k); 

old_likelihood_pre_thermo=likelihood_gen(10.^position,data_file, traj_option);
old_likelihood=thermo_temp*old_likelihood_pre_thermo;

init_likelihood=old_likelihood;
init_prior=old_prior;

old_post=old_likelihood+old_prior;
time=num2str(ceil(20.*rand));

 output_root='/files/ImStor/sorger/data/computation/Hoda';
 output_name=[output_root '/' output];

%output_name=output;

for p=start_p:time_max
    if mod(p,10000)==0
        disp(p)
        disp(output_name)
    end
    if mod(p,10000)==0
        save(output_name);     
    end
    switch lower(hess_option)
        case 'hess'
            if mod(p,first_hessian)==0
                divided=p/first_hessian;
                if p~=time_max
                    %[value index] = min(old_posterior(1,p-999:p-1));
                    %position=(new_pos(index,1:78));
                    switch traj_option
                        case 1
                            disp('case 1, hess_parallel_nodes call')
                            [final_hess_c]=hess_parallel_nodes(position,num_nodes,traj_option,prior_option,thermo_temp,k,data_file);
                        case 2
                            disp('case 2, hess_parallel_nodes call')
                            data=[data_ec,data_ic];
                            [final_hess_c]=hess_parallel_nodes(position,num_nodes,traj_option,prior_option,thermo_temp,k,data_file);
                        case 3
                            disp('case 3, hess_parallel_nodes call')
                            data=[data_ec,data_ic];
                            [final_hess_c]=hess_parallel_nodes(position,num_nodes,traj_option,prior_option,thermo_temp,k,data_file);
                    end
                    all_hessians{divided+1}=final_hess_c;
                    
                    fprintf('-----------HESSIAN CALCULATED------------ \n');
                end
            end
    end



    switch lower(boundary_option)
        case 'boundary'
                %-------------------------BOUNDARY CHECK---------------------------
                boundary=1;
                rej_count=-1;
             while boundary>0;  

                 switch lower(hess_option)
                     case 'hess'

                        if p<first_hessian
                          new_position=rand_gauss(position, sig_value, size_k, norm_step_size,k);
                        end

                       if p>=first_hessian
                           [new_position hess_ss]=step_generator(final_hess_c, position,num_param,k);
                           hess_ss_all(p)=hess_ss;
                       end

                     case 'nohess'
                         new_position=rand_gauss(position, sig_value, size_k, norm_step_size,k);
                     
                     case 'single'
                         if p==1
                             new_position=rand_gauss(position,sig_value,size_k,norm_step_size,k);
                             all_par_sigma=ones(1,size(k,2));
                             parameter_changed=1;
                         else
                             [new_position,parameter_changed,all_par_sigma,acc_all_par]=...
                                 rand_gauss_single(position,k,parameter_changed,p,accept,all_par_sigma);
                             
                             total_acc_all_par(p,:)=acc_all_par;
                             total_all_par_sigma(p,:)=all_par_sigma;
                             total_parameter_changed(p)=parameter_changed;
                         end
                 end


                 parameters(p,:)= new_position(1,:);
                 boundary=boundary_check(new_position);

                 boundary_vector(p)=boundary;
                  rej_count=rej_count+1;
                  total_rej(p)=rej_count;
                  if rej_count>20
                      position=new_pos(p-2,:);
                  end
             end
        case 'noboundary'
            switch lower(hess_option)
                     case 'hess'

                        if p<first_hessian
                          new_position=rand_gauss(position, sig_value, size_k, norm_step_size,k);
                        end

                        if p>=first_hessian
                           [new_position hess_ss]=step_generator(final_hess_c, position,num_param,k);
                           hess_ss_all(p)=hess_ss;
                      	end
                       
                     case 'nohess'                        
                         new_position=rand_gauss(position, sig_value, size_k, norm_step_size,k);
                         
                     case 'single'
                         if p==1
                             new_position=rand_gauss(position,sig_value,size_k,norm_step_size,k);
                             all_par_sigma=ones(1,size_k);
                             parameter_changed=1;
                         else
                             [new_position,parameter_changed,all_par_sigma,acc_all_par]=...
                                 rand_gauss_single(position,k,parameter_changed,p,accept,all_par_sigma);
                             
                             total_acc_all_par(p,:)=acc_all_par;
                             total_all_par_sigma(p,:)=all_par_sigma;
                             total_parameter_changed(p)=parameter_changed;
                         end
            end
    end
    %------------------------------------------------------------------
    if mod(p,first_hessian/2)==0
       fprintf(fileid, '%s p=%d\n', datestr(now, 31), p);
    end
    
    new_prior=prior_gen(new_position,prior_option,k);
    new_likelihood_pre_thermo=likelihood_gen(10.^new_position,data_file, traj_option);
    new_likelihood=thermo_temp*new_likelihood_pre_thermo;
    
    new_post=new_prior+new_likelihood;
    
    delta_post=new_post-old_post;
    delta_post_matrix(p)=delta_post;
    delta_lkls(p)=new_likelihood-old_likelihood;
    delta_prior(p)=new_prior-old_prior;
    
    test_priors(p)=new_prior;

    test_lkls(p)=new_likelihood;
    
    test(p)=exp(-delta_post);
    test_position(p,:)=new_position;
    %------------------METROPOLIS-HASTINGS ALGORITHM-------------------
    
    if exp(delta_post)<1 
        old_post=new_post;
        old_likelihood_pre_thermo=new_likelihood_pre_thermo;
        old_likelihood=new_likelihood;
        old_prior=new_prior;
        position=new_position;
        acceptance=acceptance+1;
        accept(p)=1;
    else
        alpha=random('unif',0,1);
        alpha_matrix(p)=alpha;       
        if exp(-delta_post/T) > alpha
        % if exp(-delta_post)>alpha
            old_post=new_post;
            old_likelihood_pre_thermo=new_likelihood_pre_thermo;
            old_likelihood=new_likelihood;
            old_prior=new_prior;
            position=new_position;
            acceptance=acceptance+1;
            accept(p)=1;
        else
            rej_post(p)=new_post;
            accept(p)=0;
            reject(p)=1;
        end
    end
    
    
    %-----------------------------------------------------------------
    
    %-------ADJUSTING SIGMA & TEMPERATURE (ANNEALING)--------%
    if rem(p,25)==0 %2500
        if acceptance/p <0.3
            if sig_value>0.25
                sig_value=sig_value-0.25/2;
            end
        else 
            if sig_value<=1
            sig_value=sig_value+0.25/2;
            end
        end
    end
    sigma_matrix(p)=sig_value;
    
    T=T_init* exp(-p/(first_hessian/6));%7500);
    if or(p>=first_hessian,strcmp(simulated_annealing,'no'));
        T=1;
    end

    %--------------------------------------------------------%
    
    d_post(p)=delta_post; 
    old_posterior(p)=old_post;
    new_pos(p,:)=position;
  
    
    priors(p)=old_prior;
    lkls(p)=old_likelihood; 
    lkls_pre_thermo(p)=old_likelihood_pre_thermo;
    total_time(p)=toc;
    
%     if mod(p,100)==0
%         [value index] = max(old_posterior);
%         opt_1=(new_pos(index-1,1:78));
%         figure
%         [actual_v actual_kdeg actual_kc, actual_kf, actual_kr] = values;
%         actual_pos(1:31)=actual_kf;
%         actual_pos(32:62)=actual_kr; 
%         [t_opt, y_opt]=jacobian(0,data_t,10.^new_pos(index-1,1:78));
%         [t_init, y_init]=jacobian(0,data_t,10.^init_pos);
%         plot(data_t,(data),'-r');
%         hold on;
%         plot(data_t,(y_init(:,23)/1e6),'b');
%         plot(data_t,(y_opt(:,23)/1e6),'k');
%         legend('actual pos', 'init pos', 'opt 1')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
%     end
        
end
toc

s=size(new_pos,1);
np1=new_pos(1:s/4,:);
np2=new_pos(s/4+1:s/2,:);
np3=new_pos(s/2+1:3*s/4,:);
np4=new_pos(3*s/4+1:s,:);

save(output_name);
fclose(fileid);


% figure;
% plot(1:p,old_posterior);
% figure;
% [value index]=min(old_posterior);


% % %---------------------DATA VISUALIZATION--------------------%
% %optimal value at old_obj_fn(n) & new_pos(n-1)
% [actual_v actual_kdeg actual_kc, actual_kf, actual_kr] = values;
% actual_pos(1:31)=actual_kf;
% actual_pos(32:62)=actual_kr;    
% n=63;
% for m=[1,3,5:10,12,19:21,23,25,29,30]
%     actual_pos(n)=actual_kc(m);
%     n=n+1;
% end;
% figure
% old_posterior=old_posterior';
% plot((old_posterior));
% title('objective function trajectory for the 78 parameters')
% xlabel('monte carlo step number')
% ylabel('-log_1_0 (posterior value)')
% 
% % figure;
% % plot(exp(old_posterior));
% 
% [value index] = min(old_posterior);
% opt_1=(new_pos(index,1:78));
% figure
% plot(log10(actual_pos),'*-r');
% hold on;
% plot(init_pos,'*-b')
% plot(opt_1,'*-k');
% legend('actual pos', 'init pos', 'opt 1')
% 
% switch traj_option
%     case 1
%         [t_opt, y_opt]=jacobian(0,data_t,10.^new_pos(index-1,1:78));
%         [t_init, y_init]=jacobian(0,data_t,10.^init_pos);
%         figure
%         plot(data_t(1:112,:),(data(1:112,:)),'or');
%         hold on;
%         plot(data_t(1:112,:),(y_init(1:112,23)/1e6),'b');
%         plot(data_t(1:112,:),(y_opt(1:112,23)/1e6),'ok');
%         legend('actual position', 'initial position', 'optimal position')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
%     case 2
%         [t_opt, y_opt]=jacobian(0,data_t,10.^new_pos(index,1:78));
%         [t_init, y_init]=jacobian(0,data_t,10.^init_pos);
%         figure
%         plot(data_t,(data_ec),'or');
%         hold on;
%         plot(data_t,(y_init(:,23)/1e6),'b');
%         plot(data_t,(y_opt(:,23)/1e6),'k');
%         legend('actual pos', 'init pos', 'opt 1')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
% 
%         figure;
%         plot(data_t,(data_ic),'or');
%         hold on;
%         plot(data_t,((y_init(:,26)+y_init(:,28)+y_init(:,30))/6e4),'b');
%         plot(data_t,((y_opt(:,26)+y_opt(:,28)+y_opt(:,30))/6e4),'k');
%         legend('actual pos', 'init pos', 'opt 1')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
%     case 3
%         [t_opt, y_opt]=jacobian(0,data_t,10.^new_pos(index,1:78));
%         [t_init, y_init]=jacobian(0,data_t,10.^init_pos);
%         figure
%         plot(data_t,(data_ec),'or');
%         hold on;
%         plot(data_t,(y_init(:,23)/1e6),'b');
%         plot(data_t,(y_opt(:,23)/1e6),'k');
%         legend('actual pos', 'init pos', 'opt 1')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
% 
%         figure;
%         plot(data_t,(data_ic),'or');
%         hold on;
%         plot(data_t,((y_init(:,26)+y_init(:,28)+y_init(:,30))/6e4),'b');
%         plot(data_t,((y_opt(:,26)+y_opt(:,28)+y_opt(:,30))/6e4),'k');
%         legend('actual pos', 'init pos', 'opt 1')
%         title('cPARP trajectory');
%         xlabel('time')
%         ylabel('cPARP concentration');
%         
%         figure;
%         plot(data_t, data_ec,'or');
%         hold on;
%         plot(data_t,data_ic,'og');
%         plot(data_t, (y_init(:,55)+y_init(:,57))/9e4,'b');
%         plot(data_t, (y_opt(:,55)+y_opt(:,57))/9e4,'k');

end
%----------------------------------------------------------------%