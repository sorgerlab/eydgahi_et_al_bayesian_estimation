 function lkl_for_all_kd=likelihood_gen(position,data_file, traj_option)
 load(data_file)
 %data_t=data_t*60; %david andrews data is in min, not seconds.
 %data=data(:,2:end); % we dont care about liposome data;
 
    global model_ode_observables init_cond_file
    data_existance= exist('data');
    if data_existance==0
        data=data_ec;
        num_traj=traj_option;
    end
    num_kd=size(data,2);
    %data_t=data_file(:,1);
    [ode_observables, kd_values, kd_index ic_index]=model_ode_observables();
    conc=init_cond_file();
    
    parfor_index_matrix=zeros(2, traj_option*size(kd_values,1));
    for i=1:traj_option
        parfor_index_matrix(1,(i-1)*size(kd_values,1)+1:i*size(kd_values))=i;
        parfor_index_matrix(2,(i-1)*size(kd_values,1)+1:i*size(kd_values))=1:size(kd_values,1);
    end
    if ~isempty(kd_values)
        for i=1:size(parfor_index_matrix,2)
            traj=parfor_index_matrix(1,i);
            kd=parfor_index_matrix(2,i);
            
            position_kd=position;
            position_kd(kd_index{kd})=kd_values{kd};
            kd_data=data(:,kd);
            kd_conc=conc;
            kd_conc(ic_index{kd})=kd_values{kd};
            sigma_data=0.25*kd_data;
            sigma_data((sigma_data==0))=0.25;
            [model_t, model_parts]=jacobian([],data_t,position_kd,[],[],kd_conc);
            model=sum(model_parts(:,ode_observables{traj}),2);

            likelihood_diff(i,:)=1/2*((model-kd_data).^2)./(sigma_data.^2);
            sum_lkl(i)=sum(likelihood_diff(kd,:),2);
            likelihood(i)=size(kd_data,1)/2*log(2*pi)+sum(log(sigma_data))+sum(likelihood_diff(kd,:),2);

               plot(data_t,kd_data,'-ob');hold on; plot(data_t,model,'-r');
        end
                %disp(likelihood_diff)
                %disp(sum_lkl)
                %disp(likelihood)
        lkl_for_all_kd=sum(likelihood);
                %clear likelihood;
        
    elseif isempty(kd_values)
        [model_t, model_parts]=jacobian([],data_t,position,[],[],conc);
        
        if size(model_parts,1)~=241
            save('ode_solve_error.mat')
        end
        
        model_ec=sum(model_parts(:,ode_observables{1}),2);
        model_ec=model_ec./1e6;

        remainder=mod(momp_start,3);
        r=1;

        for q=1:3:momp_start-remainder-1
            model_short_ec(r,1)=model_ec(q,1);
            data_short_ec(r,1)=data_ec(q,1);       
            time_short(r,1)=data_t(q,1);
            sigma_short_ec(r,1)=sigma_ec(q,1);
            r=r+1;
        end

        model_short_ec(r:momp_end-momp_start+r)=model_ec(momp_start:momp_end,1); 
        data_short_ec(r:momp_end-momp_start+r)=data_ec(momp_start:momp_end,1);    
        time_short(r:momp_end-momp_start+r)=data_t(momp_start:momp_end);
        sigma_short_ec(r:momp_end-momp_start+r)=sigma_ec(momp_start:momp_end);

        likelihood_diff_ec=1/2*((model_short_ec-data_short_ec).^2); %finds only the objective function for cPARP
       % likelihood_ec=size(data_short_ec,1)/2*log(2*pi)+sum(log(sigma_short_ec))+sum(likelihood_diff_ec./(sigma_short_ec.^2)); %finds the ln[p(data|k)] or the ln(likelihood)
        likelihood_ec=sum(likelihood_diff_ec./(sigma_short_ec.^2));
%                 likelihood_diff(traj,:)=1/2*((model-data).^2)./(sigma_data.^2);
%                 sum_lkl(traj)=sum(likelihood_diff(traj,:),2);
%                 likelihood(traj)=size(data,1)/2*log(2*pi)+sum(log(sigma_data))+sum(likelihood_diff(traj,:),2);
  %            plot(data_t,data,'-ob');hold on; plot(data_t,model,'-r');

%                 disp(likelihood_diff)
%                 disp(sum_lkl)
%                 disp(likelihood)
%                 lkl_for_all_kd(traj)=sum(likelihood);
                %clear likelihood;

        if or(traj_option==2, traj_option==3)
            data_ic=data(:,2);

            sigma_ic(1:momp_start-1)=0.0482;%0.0034;
            sigma_ic(momp_start:momp_end)=0.1102;%0.0227;
            sigma_ic(momp_end+1:size(data_ic))=0.0482;%0.0034;
            sigma_ic=sigma_ic';

            model_ic(:,1)=sum(model_parts(:,ode_icrp));
            model_ic=model_ic./6e4;

            r=1;
            for q=1:3:momp_start-remainder+1
                model_short_ic(r,1)=model_ic(q,1);
                data_short_ic(r,1)=data_ic(q,1);
                sigma_short_ic(r,1)=sigma_ic(q,1);

                r=r+1;
            end

            model_short_ic(r:momp_end-momp_start+r)=model_ic(momp_start:momp_end,1);
            data_short_ic(r:momp_end-momp_start+r)=data_ic(momp_start:momp_end,1);
            sigma_short_ic(r:momp_end-momp_start+r)=sigma_ic(momp_start:momp_end);

            likelihood_diff_ic=1/2*((model_short_ic-data_short_ic).^2); %finds only the objective function for cPARP
            likelihood_ic=sum(likelihood_diff_ic./(sigma_short_ic.^2)); %finds the ln[p(data|k)] or the ln(likelihood)
        end

        if num_traj==1
            likelihood=likelihood_ec;
        elseif num_traj==2
            likelihood=likelihood_ec+likelihood_ic;
        end
        lkl_for_all_kd=likelihood;
    end
   
    
%         %FOR PAPER RESUBMISSION FIGURES ONLY!*
%      model_short_ec(28:41,1)=model_ec(60:3:99,1);
%     data_short_ec(28:41,1)=data_ec(60:3:99,1);       
%      time_short(28:41,1)=data_t(60:3:99,1);
%     sigma_short_ec(28:41,1)=sigma_ec(60:3:99,1);
%     %END%

%     figure; plot(data_t,data,'ob'); hold on; plot(data_t,model_ec,'-r')
%     figure; 
%     plot(time_short,data_short_ec,'og'); hold on; plot(time_short,model_short_ec,'-b')
 end
