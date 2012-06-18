 function position=albeck_nominal_values()

    v = 0.01;
 
    kf(1)    = 4E-7;
    kf(2)    = 1E-6;
    kf(3)    = 1E-7;  %1E-6;
    kf(4)    = 1E-6;
    kf(5)    = 1E-7;
    kf(6)    = 1E-7;  %5E-7; %1E-6;
    kf(7)    = 1E-7;  %3E-8;
    kf(8)    = 2E-6;
    kf(9)    = 1E-6;
    kf(10)   = 1E-7;
    kf(11)   = 1E-6;
    kf(12)   = 1E-7;
    kf(13)   = 0.01;  % p_other(transloc);
    kf(14)   = 1E-6;
    kf(15)   = 1E-6;
    kf(16)   = 1E-6;
    kf(17)   = 1E-6;
    kf(18)   = 1E-6;
    kf(19)   = 1E-6;
    kf(20)   = 2E-6;
    kf(21)   = 2E-6;
    kf(22)   = 0.01/v;  % p_other(transloc)/p_other(v);
    kf(23)   = 5E-7;
    kf(24)   = 5E-8;
    kf(25)   = 5E-9;
    kf(26)   = 0.01/v;  % p_other(transloc)/p_other(v);
    kf(27)   = 2E-6;
    kf(28)   = 7E-6;
    kf(29)	 = 5E-8;
    kf(30)   = 6E-7;
    kf(31)   = 1E-2;  %1E-3;
    
% REVERSE RATES
    kr(1)    = 1E-6;  %1E-3;
    kr(2)    = 1E-3;
    kr(3)    = 1E-3;
    kr(4)    = 1E-3;
    kr(5)    = 1E-3;
    kr(6)    = 1E-3;
    kr(7)    = 1E-3;
    kr(8)    = 1E-3;
    kr(9)    = 1E-3;
    kr(10)   = 1E-3;
    kr(11)   = 1E-3;
    kr(12)   = 1E-3;
    kr(13)   = 0.01/v;  % p_other(transloc)/p_other(v);
    kr(14)   = 1E-3;
    kr(15)   = 1E-3;
    kr(16)   = 1E-3;
    kr(17)   = 1E-3;
    kr(18)   = 1E-3;
    kr(19)   = 1E-3;
    kr(20)   = 1E-3;
    kr(21)   = 1E-3;
    kr(22)   = 0.01;  % p_other(transloc);
    kr(23)   = 1E-3;
    kr(24)   = 1E-3;
    kr(25)   = 1E-3;
    kr(26)   = 0.01;  % p_other(transloc);
    kr(27)   = 1E-3;
    kr(28)   = 1E-3;
    kr(29)	 = 1E-3;
    kr(30)   = 1E-3;
    kr(31)   = 1e-3;    %0;    
    
% "catalysis" rates
    kc(1)    = 1E-2;  %5E-5; % 1E-5;
    kc(3)    = 1;   
    kc(5)    = 1;   
    kc(6)    = 1;   
    kc(7)    = 1;   
    kc(8)    = 0.1;  %CC3UBRate; 
    kc(9)    = 20;  %1;     
    kc(10)   = 1;
    kc(12)   = 1;
    kc(19)   = 1;
    kc(20)   = 10;  
    kc(21)   = 10;  
    kc(23)   = 1;   
    kc(25)   = 1;
    kc(29)   = 1;
    kc(30)   = 1e-3;    %0;  
    
    position(1:31)=kf;
    position(32:62)=kr;
    
    n=63;
    for m=[1,3,5:10,12,19:21,23,25,29,30]
        position(n)=kc(m);
        n=n+1;
    end;
    
%     if size(mult_factor,2)==96
%     
%         conc = values_conc;
% 
%         %have covariance data for the following 5 IC's
%         position(79)=conc(33); %Bcl2
%         position(80)=conc(24); %Bid
%         position(81)=conc(12); %C3
%         position(82)=conc(19); %PARP
%         position(83)=conc(29); %Bax
% 
%         position(84)=conc(5);
%         position(85)=conc(7);
%         position(86)=conc(10);
%         position(87)=conc(15);
%         position(88)=conc(21);
%         position(89)=conc(27);
%         position(90)=conc(39);
%         position(91)=conc(42);
%         position(92)=conc(45);
%         position(93)=conc(49);
%         position(94)=conc(52);
% 
%         %ligand & receptor
%         position(95)=conc(2); %receptor
%         position(96)=conc(1); %ligand
%     end
  
   % position=position.*mult_factor;
 end