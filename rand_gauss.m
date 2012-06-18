function log_new_position=rand_gauss(position, sig_value, size_k, norm_step_size,k)  
%generates the log10 of the random new positions

%this m-file uses a random gaussian to generate a random step size (centered 
%at the origin with a sigma of 1) and adds it to the previous position in 
%log space to generate the new position of the parameters

%this m-file is used by MH.m
    mu=zeros(1,size_k);%log(k)
    sigma=sig_value*ones(1,size_k);
    obj = gmdistribution(mu, sigma);
    step_position = random(obj);
    length=norm(step_position);
    delta_position=step_position./length;
    if size_k==2
        delta_position(1,3:78)=0;
    elseif size_k==96
        delta_position(1,end)=0;
    elseif size(position,2)~=size_k
        %num_ic=size(position,2)-size_k;
        delta_position_temp=zeros(1,size(position,2));
        delta_position_temp(1,k)=delta_position;
        clear delta_position;
        delta_position=delta_position_temp;
    end
        
    log_new_position=position+norm_step_size*delta_position;
end