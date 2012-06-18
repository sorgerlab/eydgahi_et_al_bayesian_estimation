function [new_step hess_ss]=step_generator(hessian, current_step,num_param,k)

%generates a new step when exploring parameter space in the direction of
%eigenvalues/eigenvectors
[eig_vec, eig_val]=eig(hessian);
num_k=size(hessian,2);

hess_ss=0.085;

ev=diag(eig_val);
for i=1:num_k
    if abs(ev(i))<0.25
        ev(i)=0.25;
    end
end
total_step=zeros(1,num_k);
for i=1:num_k
    sigma=1./sqrt(abs(ev(i)));
    mu=0;%zeros(1,s);%log(k)
    obj = gmdistribution(mu, sigma);
    step_position = random(obj);
    delta_step=step_position*hess_ss*eig_vec(:,i)'; %*0.01
    
    total_step=total_step+delta_step;
end

%     if num_parm==96
%         total_step(1,end)=0;
%     elseif num_parm==2
%         total_step(1,3:size(current_step,2))=0;
    if size(current_step,2)~=size(num_k,2)
        %num_ic=size(current_step,2)-num_param;
        %total_step_temp(1,1:num_ic)=0;
        %total_step_temp(1,num_ic+1:num_ic+size(current_stop,2))=total_step;
        %clear total_step;
        
        total_step_temp=zeros(1,num_param);
        total_step_temp(k)=total_step;
        clear total_step;
        total_step=total_step_temp;
    end
    
new_step=current_step+total_step;
end