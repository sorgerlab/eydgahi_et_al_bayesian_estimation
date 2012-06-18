function [prior]=prior_gen(position, prior_option, k)

global prior_flag_file

prior_flag=prior_flag_file();
%not a rate contant: 0
%forward, 1st order: 1
%reverse, 1st order: 2
%foward, 2nd order: 3
%reverse, 2nd order: 4
%catalytic: 5

s=size(prior_flag,2);
prior=zeros(1,s);

switch lower(prior_option)
    case 'uniform'
        prior=0;
    case 'albeck'
        [u_f1,u_f2,u_r1,u_r2,u_c,v_f1,v_f2,v_r1,v_r2,v_c]=prior_stats;
        for i=1:s
            if prior_flag(i)==1
                prior(i)=(1/sqrt(2*pi*v_f1))*exp(-(position(i)-u_f1)^2/(2*v_f1));
            elseif prior_flag(i)==3
                prior(i)=(1/sqrt(2*pi*v_f2))*exp(-(position(i)-u_f2)^2/(2*v_f2));
            elseif prior_flag(i)==2
                prior(i)=(1/sqrt(2*pi*v_r1))*exp(-(position(i)-u_r1)^2/(2*v_r1));
            elseif prior_flag(i)==4
                prior(i)=(1/sqrt(2*pi*v_r2))*exp(-(position(i)-u_r2)^2/(2*v_r2));
            elseif prior_flag(i)==5
                prior(i)=(1/sqrt(2*pi*v_c))*exp(-(position(i)-u_c)^2/(2*v_c));
            elseif prior_flag(i)==0
                prior(i)=0;
            end
        end
        prior=-1*sum(log(prior(k))); %sum over the ln(prior(k_i)). multiplying 78 numbers is very small. ln of it... not so bad.
   
        
    case 'general'
        [u_f1,u_f2,u_r1,u_r2,u_c,v_f1,v_f2,v_r1,v_r2,v_c]=prior_stats_old;       
        for i=1:s
            if prior_flag(i)==1
                prior(i)=(1/sqrt(2*pi*v_f1))*exp(-(position(i)-u_f1)^2/(2*v_f1));
            elseif prior_flag(i)==3
                prior(i)=(1/sqrt(2*pi*v_f2))*exp(-(position(i)-u_f2)^2/(2*v_f2));
            elseif prior_flag(i)==2
                prior(i)=(1/sqrt(2*pi*v_r1))*exp(-(position(i)-u_r1)^2/(2*v_r1));
            elseif prior_flag(i)==4
                prior(i)=(1/sqrt(2*pi*v_r2))*exp(-(position(i)-u_r2)^2/(2*v_r2));
            elseif prior_flag(i)==5
                prior(i)=(1/sqrt(2*pi*v_c))*exp(-(position(i)-u_c)^2/(2*v_c));
            elseif prior_flag(i)==0
                prior(i)=0;
            end
        end
        prior=-1*sum(log(prior(k))); %sum over the ln(prior(k_i)). multiplying 78 numbers is very small. ln of it... not so bad.
end

end

