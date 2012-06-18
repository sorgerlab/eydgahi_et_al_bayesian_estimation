function [oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,element_num,...
    position_ff_total, position_fb_total,position_bf_total,position_bb_total]=...
    hess_gradient_calc(position,ii,jj,diagonal,element_num,s,...
    oneD_hessian,back_hessian,hess_ff,hess_fb,hess_bf,hess_bb,...
    traj_option,prior_option,thermo_temp,k, data_file)    


total_elements=sum(1:s);
delta_k=0.1;

position_ff=position;
position_fb=position;
position_bf=position;
position_bb=position;

for i=ii
    position_ff(k(i))=position_ff(k(i))+delta_k;
    position_fb(k(i))=position_fb(k(i))+delta_k;
    position_bf(k(i))=position_bf(k(i))-delta_k;
    position_bb(k(i))=position_bb(k(i))-delta_k;
 

    for j=jj%count_2:s;
        if j==diagonal
            oneD_hessian_lkl(j)=thermo_temp*likelihood_gen(10.^position_ff,data_file, traj_option);
            back_hessian_lkl(j)=thermo_temp*likelihood_gen(10.^position_bb,data_file, traj_option);

            oneD_hessian_prior(j)=prior_gen(position_ff,prior_option, k);
            back_hessian_prior(j)=prior_gen(position_bb,prior_option, k);

            oneD_hessian(j)=oneD_hessian_lkl(j)+oneD_hessian_prior(j);
            back_hessian(j)=back_hessian_lkl(j)+back_hessian_prior(j);
            
            
            position_ff_total{element_num,1}=position_ff;
            position_fb_total{element_num,1}=position_fb;
            position_bf_total{element_num,1}=position_bf;
            position_bb_total{element_num,1}=position_bb;
            element_num=element_num+1;
            
            [i,j, oneD_hessian(j),element_num-1]
        else
            
            position_ff(k(j))=position_ff(k(j))+delta_k;
            position_fb(k(j))=position_fb(k(j))-delta_k;
            position_bf(k(j))=position_bf(k(j))+delta_k;
            position_bb(k(j))=position_bb(k(j))-delta_k;
            
            hess_ff_lkl(i,j)=thermo_temp*likelihood_gen(10.^position_ff,data_file, traj_option);
            hess_fb_lkl(i,j)=thermo_temp*likelihood_gen(10.^position_fb,data_file, traj_option);
            hess_bf_lkl(i,j)=thermo_temp*likelihood_gen(10.^position_bf,data_file, traj_option);
            hess_bb_lkl(i,j)=thermo_temp*likelihood_gen(10.^position_bb,data_file, traj_option);

            hess_ff_prior(i,j)=prior_gen(position_ff,prior_option,k);
            hess_fb_prior(i,j)=prior_gen(position_fb,prior_option,k);
            hess_bf_prior(i,j)=prior_gen(position_bf,prior_option,k);
            hess_bb_prior(i,j)=prior_gen(position_bb,prior_option,k);

            hess_ff(element_num)=hess_ff_lkl(i,j)+hess_ff_prior(i,j);
            hess_fb(element_num)=hess_fb_lkl(i,j)+hess_fb_prior(i,j);
            hess_bf(element_num)=hess_bf_lkl(i,j)+hess_bf_prior(i,j);
            hess_bb(element_num)=hess_bb_lkl(i,j)+hess_bb_prior(i,j);
            
            position_ff_total{element_num,1}=position_ff;
            position_fb_total{element_num,1}=position_fb;
            position_bf_total{element_num,1}=position_bf;
            position_bb_total{element_num,1}=position_bb;

            position_ff(k(j))=position_ff(k(j))-delta_k;
            position_fb(k(j))=position_fb(k(j))+delta_k;
            position_bf(k(j))=position_bf(k(j))-delta_k;
            position_bb(k(j))=position_bb(k(j))+delta_k;
           
            
            [i,j, hess_ff(element_num), hess_fb(element_num), element_num]
            element_num=element_num+1;
            
        end
    end
    position_ff(k(i))=position_ff(k(i))-delta_k;
    position_fb(k(i))=position_fb(k(i))-delta_k;
    position_bf(k(i))=position_bf(k(i))+delta_k;
    position_bb(k(i))=position_bb(k(i))+delta_k;
end
end