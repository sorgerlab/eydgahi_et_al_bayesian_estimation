% function output=hessian(position)
%-------------initialization------------------
function [final_hess_c,t,oneD_hessian, back_hessian, hess_ff, hess_fb, hess_bf, hess_bb, orig_value,...
    position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
    hessian_parallel_short(start,finish,job_number, position,row_number,traj_option,prior_option,thermo_temp,k, model_ode_observables, ode_file, prior_flag_file, init_cond_file, data_file, likelihood_gen_file)
tic
global model_ode_observables 
global ode_file 
global prior_flag_file 
global init_cond_file
global likelihood_gen_file

s=size(k,2); 

orig_prior=prior_gen(position,prior_option,k);
orig_lkl=thermo_temp*likelihood_gen(10.^position,data_file, traj_option);
orig_value=orig_lkl+orig_prior;
%s=size(position); s=s(2);

clear final_hess_c, clear actual_delta,
clear hess_*, clear gradient_c, clear predicted_delta_*

total_elements=sum(1:s);
delta_k=0.1;

oneD_hessian=zeros(1,s);                       
back_hessian=zeros(1,s);
hess_ff=zeros(1,total_elements);
hess_fb=zeros(1,total_elements);
hess_bf=zeros(1,total_elements);
hess_bb=zeros(1,total_elements);

start_row=find(start>=row_number,1,'last');
finish_row=find(finish>=row_number,1,'last');

start_column=mod(start,row_number(start_row))+start_row;
finish_column=mod(finish,row_number(finish_row))+finish_row;

if start_row==1
    start_column=start;
end

if finish_row==1
    finish_column=finish;
    start_column=start;
end

element_num=start;
if start_row==finish_row
    i=start_row;
    j=start_column:finish_column;
    [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
        position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
        hess_gradient_calc(position,i,j,start_row,element_num,s,...
        oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
        traj_option,prior_option,thermo_temp,k, data_file); 
else
    if start_row==finish_row-1
        i=start_row:finish_row-1;
        j=start_column:s;
        [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
            position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
            hess_gradient_calc(position,i,j,i,element_num,s,...
            oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
            traj_option,prior_option,thermo_temp,k, data_file); 
    else
        i=start_row;
        j=start_column:s;
        [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
            position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
            hess_gradient_calc(position,i,j,i,element_num,s,...
            oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
            traj_option,prior_option,thermo_temp,k, data_file);
        
        for i=start_row+1:finish_row-1;
            j=i:s;
            [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
                position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
                hess_gradient_calc(position,i,j,i,element_num,s,...
                oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
                traj_option,prior_option,thermo_temp,k, data_file);
        end

    end
    i=finish_row;
    j=finish_row:finish_column;
    [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
        position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
        hess_gradient_calc(position,i,j,i,element_num,s,...
        oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
        traj_option,prior_option,thermo_temp,k, data_file);
    
      
end
    
    for element_num=start:finish
        if ismember(element_num,row_number)
            j=find(element_num==row_number);
            final_hess_c(element_num)=(oneD_hessian(j)-2*orig_value+back_hessian(j))/delta_k^2;
        else
            final_hess_c(element_num)=(hess_ff(element_num)-hess_fb(element_num)-hess_bf(element_num)+hess_bb(element_num))/(4*delta_k^2);
        end
    end
        
    toc
    t=toc;
    save_name=['hess_par_',num2str(job_number),'.mat'];
    save(save_name)
end