function conc=albeck_init_conds()
global v kdeg ks
    conc0_c(2)    = 3000;%3E3;  % 1000ng/ml=60000molec/cell, 250 = 15000 50 = 3000, 10 = 600  % Ligand
    conc0_c(1)    = 1E3;  % WT = 1E3; RNAi = 2E3;   % Receptor
    conc0_c(35)   = 2E3;  % John's model has 1E2;  used 2E3 so 1E2 = floor with 20x up/down variation  % flip
    conc0_c(27)   = 1E4;  % WT = 1E4; RNAi = 2E4;   % C8
    conc0_c(320)  = 1E3;  % Bar
    conc0_c(57)   = 1E4;  % C3
    conc0_c(322)  = 1E4;  % C6
    conc0_c(108)  = 1E5;  % XIAP
    conc0_c(105)  = 1E6;  % PARP
    conc0_c(181)  = 6E4;  % WT = 6E4; RNAi = 1.2E5; % Bid
    conc0_c(326)  = 2E4;  % WT = 2E4; RNAi = 2E4;   % Mcl1 (Bclx)
    conc0_c(182)  = 8E4;  % WT = 8E4; RNA1 = 4E5;   % Bax
    conc0_c(206)  = 3E4;  %Feb09 change 6E4;  % Bcl2
    conc0_c(334)  = 5E5;  % M
    conc0_c(188)  = 5E5;  % CyCm
    conc0_c(339)  = 1E5;  % Smac_m
    conc0_c(191)  = 1E5;  % Apaf
    conc0_c(193)  = 1E5;  % C9
    conc0_c(348)  = 0;    %IETD
    conc0_c(351)  = 0;    %DEVDR
    
 %initial values
    conc(1)     = conc0_c(2);
    conc(2)     = conc0_c(1);
    conc(3:4)   = 0;
    conc(5)     = conc0_c(35);
    conc(6)     = 0;
    conc(7)     = conc0_c(27);
    conc(8:9)   = 0;
    conc(10)    = conc0_c(320);
    conc(11)    = 0;
    conc(12)    = conc0_c(57);
    conc(13:14) = 0;
    conc(15)    = conc0_c(322);
    conc(16:18) = 0;
    conc(19)    = conc0_c(108);
    conc(20)    = 0;
    conc(21)    = conc0_c(105);
    conc(22:23) = 0;
    conc(24)    = conc0_c(181);
    conc(25:26) = 0;
    conc(27)    = conc0_c(326);
    conc(28)    = 0;
    conc(29)    = conc0_c(182);
    conc(30:32) = 0;
    conc(33)    = conc0_c(206);
    conc(34:38) = 0;
    conc(39)    = conc0_c(334);
    conc(40:41) = 0;
    conc(42)    = conc0_c(188);
    conc(43:44) = 0;
    conc(45)    = conc0_c(339);
    conc(46:48) = 0;
    conc(49)    = conc0_c(191);
    conc(50:51) = 0;
    conc(52)    = conc0_c(193);
    conc(53:69) = 0;
    
    
    v = 0.01;
     kdeg_id(2)   = 2.9E-6; %1/2*1/24*1/3600;  
    kdeg_id(1)   = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(3)   = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(317) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(35)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(318) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(27)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(319) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(56)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(320) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(321) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(57)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(58)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(59)  = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(322) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(323) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(324) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(325) = 2.9E-6; %1/2*1/24*1/3600;     
    kdeg_id(108) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(210) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(105) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(106) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(107) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(181) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(183) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(184) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(326) = 1E-4; %1/2*1/24*1/3600;
    kdeg_id(327) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(182) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(185) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(328) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(329) = 2.9E-6; %1/2*1/24*1/3600;     
    kdeg_id(206) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(330) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(187) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(331) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(332) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(333) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(334) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(335) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(336) = 1E-4;
    kdeg_id(188) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(337) = 2.9E-6; %1/2*1/24*1/3600;  
    kdeg_id(338) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(339) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(340) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(341) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(190) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(191) = 2.9E-6; %1/2*1/24*1/3600; 
    kdeg_id(192) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(342) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(193) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(343) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(344) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(345) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(346) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(347) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(223) = 0;    
    kdeg_id(348) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(349) = 2.9E-6; %1/2*1/24*1/3600;     
    kdeg_id(350) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(351) = 2.9E-6; %1/2*1/24*1/3600;    
    kdeg_id(352) = 2.9E-6; %1/2*1/24*1/3600;
    kdeg_id(353) = 2.9E-6; %1/2*1/24*1/3600; 
        
    kdeg_identity=[ 'kdeg_id(2)  '; 'kdeg_id(1)  '; 'kdeg_id(3)  '; 'kdeg_id(317)'; 
        'kdeg_id(35) '; 'kdeg_id(318)'; 'kdeg_id(27) '; 'kdeg_id(319)'; 'kdeg_id(56) ';
        'kdeg_id(320)'; 'kdeg_id(321)'; 'kdeg_id(57) '; 'kdeg_id(58) '; 'kdeg_id(59) '; 
        'kdeg_id(322)'; 'kdeg_id(323)'; 'kdeg_id(324)'; 'kdeg_id(325)'; 'kdeg_id(108)'; 
        'kdeg_id(210)'; 'kdeg_id(105)'; 'kdeg_id(106)'; 'kdeg_id(107)'; 'kdeg_id(181)'; 
        'kdeg_id(183)'; 'kdeg_id(184)'; 'kdeg_id(326)'; 'kdeg_id(327)'; 'kdeg_id(182)'; 
        'kdeg_id(185)'; 'kdeg_id(328)'; 'kdeg_id(329)'; 'kdeg_id(206)'; 'kdeg_id(330)'; 
        'kdeg_id(187)'; 'kdeg_id(331)'; 'kdeg_id(332)'; 'kdeg_id(333)'; 'kdeg_id(334)'; 
        'kdeg_id(335)'; 'kdeg_id(336)'; 'kdeg_id(188)'; 'kdeg_id(337)'; 'kdeg_id(338)'; 
        'kdeg_id(339)'; 'kdeg_id(340)'; 'kdeg_id(341)'; 'kdeg_id(190)'; 'kdeg_id(191)'; 
        'kdeg_id(192)'; 'kdeg_id(342)'; 'kdeg_id(193)'; 'kdeg_id(343)'; 'kdeg_id(344)'; 
        'kdeg_id(345)'; 'kdeg_id(346)'; 'kdeg_id(347)'; 'kdeg_id(223)'; 'kdeg_id(348)'; 
        'kdeg_id(349)'; 'kdeg_id(350)'; 'kdeg_id(351)'; 'kdeg_id(352)'; 'kdeg_id(353)'];
        
    kdeg_identity=cellstr(kdeg_identity);
    
    for i=1:64
        identity_value=eval(kdeg_identity{i});
        kdeg(i)=identity_value;
    end 
    
     % synthesis rates
    ks_id(1)    = 0;    
    ks_id(2)    = 0.15*2.9E-6*conc0_c(1); %1/2*1/24*1/3600*conc0_c(1);
    ks_id(35)   = 0.15*2.9E-6*conc0_c(35); %1/2*1/24*1/3600*conc0_c(35);
    ks_id(27)   = 0.15*2.9E-6*conc0_c(27); %1/2*1/24*1/3600*conc0_c(27);
    ks_id(320)  = 0.15*2.9E-6*conc0_c(320); %1/2*1/24*1/3600*conc0_c(320); 
    ks_id(57)   = 0.15*2.9E-6*conc0_c(57); %1/2*1/24*1/3600*conc0_c(57);
    ks_id(322)  = 0.15*2.9E-6*conc0_c(322); %1/2*1/24*1/3600*conc0_c(322); 
    ks_id(108)  = 0.15*2.9E-6*conc0_c(108); %1/2*1/24*1/3600*conc0_c(108); 
    ks_id(105)  = 0.15*2.9E-6*conc0_c(105); %1/2*1/24*1/3600*conc0_c(105); 
    ks_id(181)  = 0.15*2.9E-6*conc0_c(181); %1/2*1/24*1/3600*conc0_c(181);
    ks_id(326)  = 0.15*1E-4*conc0_c(326); %1E-4*conc0_c(326);
    ks_id(182)  = 0.15*2.9E-6*conc0_c(182); %1/2*1/24*1/3600*conc0_c(182); 
    ks_id(206)  = 0.15*2.9E-6*conc0_c(206); %1/2*1/24*1/3600*conc0_c(206);
    ks_id(334)  = 0.15*2.9E-6*conc0_c(334); %1/2*1/24*1/3600*conc0_c(334);
    ks_id(188)  = 0.15*2.9E-6*conc0_c(188); %1/2*1/24*1/3600*conc0_c(188);
    ks_id(339)  = 0.15*2.9E-6*conc0_c(339); %1/2*1/24*1/3600*conc0_c(339);    
    ks_id(191)  = 0.15*2.9E-6*conc0_c(191); %1/2*1/24*1/3600*conc0_c(191); 
    ks_id(193)  = 0.15*2.9E-6*conc0_c(193); %1/2*1/24*1/3600*conc0_c(193);
	ks_id(348)  = 0.15*2.9E-6*conc0_c(348); %1/2*1/24*1/3600*conc0_c(348);
	ks_id(351)  = 0.15*2.9E-6*conc0_c(351); %1/2*1/24*1/3600*conc0_c(351);
       
       ks_identity=['ks_id(1)  '; 'ks_id(2)  '; 'ks_id(35) '; 'ks_id(27) '; 
       'ks_id(320)'; 'ks_id(57) '; 'ks_id(322)'; 'ks_id(108)';
       'ks_id(105)'; 'ks_id(181)'; 'ks_id(326)'; 'ks_id(182)';
       'ks_id(206)'; 'ks_id(334)'; 'ks_id(188)'; 'ks_id(339)'; 
       'ks_id(191)'; 'ks_id(193)'; 'ks_id(348)'; 'ks_id(351)']; 
        
    ks_identity=cellstr(ks_identity);
    for i=1:20
        index=[1,2,5,7,10,12,15,19,21,24,27,29,33,39,42,45,49,52,59,62];
        identity_value=eval(ks_identity{i});
        ks(index(i))=identity_value;
    end %this for loop assigns all the actual ks values to the indices used in the actual odes

 
end