function [t,y]=jacobian(xyz,tspan,position,kc_index,trail_conc, kd_conc)
%solves the 69 ODEs in S. Gaudet's TRAIL model and generates all 69 trajectories

%to call it jacobian.m, 
%if generating points for the actual data obtained from experiment, use: [data_t, data]=jacobian(1,[],[]);
%if generating points for the model from fake data, use: [model_t, model]=jacobian(0,data_t,new_position);

%this M-file is used by MH.m
%this M-file uses the following files:
%values_conc.m
%values.m
%odes.m
%--------------------------------------------------------------------------
global rel_tol_value nominal_value_file init_cond_file ode_file %global kc kf kr v kdeg ks rel_tol_value

if ~isempty(xyz)
    conc = init_cond_file(); % [conc, ks] = values_conc;
end
%[v, kdeg] = values;
if (xyz==1)
    tspan=0:21.6:43200; %times intervals of 21.6 used to generate 2000 evenly 
    %distributed time points for calculating objective function
    nominal_par_values= nominal_value_file();%[kc, kf, kr, v, kdeg] = values;
end

if (xyz==0)
    kf=position(1:31);
    kr=position(32:62);
    n=63;
    for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
    end;
%     if size(position,2)==96
%         conc(33)=position(79);
%         conc(24)=position(80);
%         conc(12)=position(81);
%         conc(19)=position(82);
%         conc(29)=position(83);
%         
%         conc(5)=position(84);
%         conc(7)=position(85);
%         conc(10)=position(86);
%         conc(15)=position(87);
%         conc(21)=position(88);
%         conc(27)=position(89);
%         conc(39)=position(90);
%         conc(42)=position(91);
%         conc(45)=position(92);
%         conc(49)=position(93);
%         conc(52)=position(94);
% 
%         %receptor
%         conc(2)=position(95);
        
        %dont want to include the values for ligand. want to keep that
        %constant for now
%    end


end

if (xyz==2) %for generating albeck et al, mol cell fig 5E
    kf=position(1:31);
    kr=position(32:62);
    n=63;
    for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
    end;    
    kc_8_values=[0,1e-4,3.3e-4,1e-3,3.3e-3,1e-2,1];
    
    kc(8)=kc_8_values(kc_index)*kc(8);
 %   conc(33)=30*3e4;
    
    %tf=12*3600; 
    %samp_freq=tf/60; 
    %dt=(tf-0)/samp_freq;
    %tspan=0:dt:tf;
    %options=odeset('AbsTol',1E-10,'RelTol',1E-8,'maxstep',dt);
    
end

if (xyz==3) %for generating figure for t_d & t_s (changes the value of trail)
    kf=position(1:31);
    kr=position(32:62);
    n=63;
    for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
    end;
 %   conc(1)=trail_conc;
end

if (xyz==10) %FOR KNOCKDOWNS, IT BEGINS AT 10! ORIGINAL CONDITIONS
    kf=position(1:31);
    kr=position(32:62);
    n=63;
    for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
    end;
end

if (xyz==11) %BCL2 OVEREXPRESSION SHOULD SUPPRESS APOPTOSIS FROM OCCURING
    kf=position(1:31);
    kr=position(32:62);
    n=63;
    for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
    end;
  %  conc(33)=30*3e4;
end

if isempty(xyz)
    conc=kd_conc;
end

warn_state = warning();
warning('off', 'MATLAB:ode15s:IntegrationTolNotMet');
warning('off', 'MATLAB:illConditionedMatrix');
options = odeset('RelTol', rel_tol_value);
%[t,y] = ode15s(ode_file,tspan, conc, options, {position});
ode_fn = @(t,y0) ode_file(t,y0,position);
[t,y]= ode15s(ode_fn, tspan, conc, options);
warning(warn_state);

if length(t)~=length(tspan)
    options = odeset('RelTol', 1e-6);
    [t,y] = ode15s(ode_fn, tspan, conc, options);
end
if length(t)~=length(tspan)
    options = odeset('RelTol', 1e-12);
    [t,y] =ode15s(ode_fn, tspan, conc, options);
end

end
% figure;
% plot(t,y(:,23),'r')
% title('cPARP trajectory');
% xlabel('time (sec)');
% ylabel('concentration');
% end

%-----S. GAUDET'S TRAIL MODEL THAT WAS ORIGINALLY CODED IN JACOBIAN------
% species=['c2  '; 'c1  '; 'c3  '; 'c317'; 'c35 '; 'c318'; 'c27 '; 'c319'; 'c56 ';
%     'c320'; 'c321'; 'c57 '; 'c58 '; 'c59 '; 'c322'; 'c323'; 'c324'; 'c325';
%     'c108'; 'c210'; 'c105'; 'c106'; 'c107'; 'c181'; 'c183'; 'c184'; 'c326';
%     'c327'; 'c182'; 'c185';'c328'; 'c329'; 'c206'; 'c330'; 'c187'; 'c331'; 
%     'c332'; 'c333'; 'c334'; 'c335'; 'c336'; 'c188'; 'c337'; 'c338'; 'c339';
%     'c340'; 'c341'; 'c190'; 'c191'; 'c192'; 'c342'; 'c193'; 'c343'; 'c344';
%     'c345'; 'c346'; 'c347'; 'c223'; 'c348'; 'c349'; 'c350'; 'c351'; 'c352';
%     'c353'; 'c354'; 'c355'; 'c356'; 'c357'; 'c358'];
% species=cellstr(species);
% for i=1:69
%     eval([species{i}, '=', num2str(i)]);
% end


% %DERIVED VARIABLES
% PARP_norm = conc(c105)/conc0_c(105) ; 
% Smac_m_norm = conc(c339)/conc0_c(339) ;
% allcPARP = conc(c107) + conc(c354) ; 
% cytoSMAC = conc(c341) + conc(c345) + conc(c347) ;
% allcytoSMAC = conc(c341) + conc(c345) + conc(c347) + conc(c355) ;
% 	
% 
% %DERVIED VARIABLES: SLOPES!
% slope_cPARP1 = dconc_dt(c107) ;
% slope_cPARP2 = dconc_dt(c356) ;
% slope_SMAC1 = dconc_dt(c357) ;
% slope_SMAC2 = dconc_dt(c358) ;


%------------MODEL Casp_Kinetics---------
%     Reactions included in the model - Note:  parameter values outdated, see below.  Reactions accurate.
%     % L(c2) + R(C1) <--> L:R(c3) --> R*(c317)
%       kf(1) := 4E-7; kr(1) := 1E-3; kc(1) := 1E-5;
%     % R*(c317) + flip(c35) <-->  R*:flip(c318)
%       kf(2) := 1E-6; kr(2) := 1E-3;
%     % R*(c317) + C8(c27) <--> R*:C8(c319) --> R*(c317) + C8*(c56)
%       kf(3) := 1E-6; kr(3) := 1E-3; kc(3) := 1;
%     % C8*(c56) + Bar(c320) <--> C8*:Bar(c321)
%       kf(4) := 1E-6; kr(4) := 1E-3;
%     % C8*(c56) + C3(c57) <--> C8*:C3(c58) --> C8*(c56) + C3*(c59)
%       kf(5) := 1E-7; kr(5) := 1E-3; kc(5) := 1;
%     % C3*(c59) + C6(c322) <--> C3*:C6(c323) --> C3*(c59) + C6?? c6*(c324)
%       kf(6) := 1E-6; kr(6) := 1E-3; kc(6):=1;
%     % C6*(c324) + C8(c27) <--> C6*:C8(c325) --> C6*(c324) + C8*(c56)
%       kf(7):=3E-8; kr(7):=1E-3; kc(7):=1;
%     % C3*(c59) + XIAP(c108) <--> C3*:XIAP(c210) --> C3*_Ub(c223) + XIAP(c108)
%       kf(8):=2E-6; kr(8):=1E-3;  kc(8):=0.1; 
%     % C3*(c59) + PARP(c105) <--> C3*:PARP(c106) --> C3*(c59) + cPARP(c107)
%       kf(9):=1E-6; kr(9):=1E-3; kc(9):=1;
%     % C8*(c56) + Bid(c181) <--> C8*:Bid(c183) --> C8*(c56) + tBid(c184)
%       kf(10):=1E-7; kr(10):=1E-3; kc(10):=1;
%     % tBid(c184) + Mcl1(c326) <-->  Bid:Mcl1(c327)
%       kf(11):=1E-6; kr(11):=1E-3; 
%     % tBid(c184) + Bax(c182) <--> tBid:Bax(c185) --> tBid(c184) + Bax*(c328)
%       kf(12):=1E-7; kr(12):=1E-3; kc(12):=1;
%     % Bax*(c328) <-->  Bax*_m(c329)
%       kf(13):=0.01; kr(13):=0.01;
%     % Bax*_m(c329) + Bcl2(c206) <-->  Bax*_m:Bcl2(c330)
%       kf(14):=1E-6; kr(14):=1E-3;
%     % Bax*_m(c329) + Bax*_m(c329) <--> Bax2(c187)
%       kf(15):=1E-6; kr(15):=1E-3;
%     % Bax2(c187) + Bcl2(c206) <--> Bax2:Bcl2(c331)
%       kf(16):=1E-6; kr(16):=1E-3; 
%     % Bax2(c187) + Bax2(c187) <--> Bax4(c332)
%       kf(17):=1E-6; kr(17):=1E-3;
%     % Bax4(c332) + Bcl2(c206) <--> Bax4:Bcl2(c333)
%       kf(18):=1E-6; kr(18):=1E-3; 
%     % Bax4(c332) + M(c334) <--> Bax4:M(c335) -->  M*(c336)
%       kf(19):=1E-6; kr(19):=1E-3; kc(19):=1;
%     % M*(c336) + CyCm(c188) <-->  M*:CyCm(c337) --> M*(c336) + CyCr(c338)
%       kf(20):=2E-6; kr(20):=1E-3; kc(20):=10;
%     % M*(c336) + Smac_m(c339) <-->  M*:Smac_m(c340) --> M*(c336) + Smac_r(c341)
%       kf(21):=2E-6; kr(21):=1E-3; kc(21):=10;
%     % CyCr(c338) <-->  CyC(c190)
%       kf(22):=0.01; kr(22):=0.01;
%     % CyC(c190) + Apaf(c191) <--> CyC:Apaf(c192) --> CyC(c190) + Apaf*(c342)
%       kf(23):=5E-7; kr(23):=1E-3; kc(23):=1;
%     % Apaf*(c342) + C9(c193) <--> Apop(c343)
%       kf(24):=5E-8; kr(24):=1E-3;
%     % Apop(c343) + C3(c57) <-->  Apop:C3(c344) --> Apop(c343) + C3*(c59)
%       kf(25):=5E-9; kr(25):=1E-3; kc(25):=1;
%     % Smac_r(c341) <-->  Smac(c345)
%       kf(26):=0.01; kr(26):=0.01;
%     % Apop(c343) + XIAP(c108) <-->  Apop:XIAP(c346)
%       kf(27):=2E-6; kr(27):=1E-3;
%     % Smac(c345) + XIAP(c108) <-->  Smac:XIAP(c347)
%       kf(28):=7E-6; kr(28):=1E-3;
%	  % IETD(c348) + C8*(c56) <--> IETD:C8*(c349) --> IETD*(c350) + C8*(c56)
%		kf(29):= 5E-8 ; kr(29):= 1e-3; kc(29):= 1;
%	  % DEVDR(c351) + C3*(c59) <--> DEVDR:C3*(c352) --> DEDVR*(c353) + C3*(c59)
%       kf(30):= 6E-7 ; kr(30):= 1E-3 ; kc(30):= 1;
%     % R*(c3) <--> L(c2) + R(c1)
%		kf(31):= 1E-3; kr(31):= 0; kc(31):= 0

%	  % L(c2) <--> null
%		ks2:=  ;  kdeg2:=  ;
%	  % R(c1) <--> null
%		ks1:=  ;  kdeg1:=  ;
%	  % flip(c35) <--> null
%		ks35:=  ;  kdeg35:=  ;
%	  % C8(c27) <--> null
%		ks27:=  ;  kdeg27:=  ;
%	  % Bar(c320) <--> null
%		ks320:=  ;  kdeg320:=  ;
%	  % C3(c57) <--> null
%		ks57:=  ;  kdeg57:=  ;
%	  % C6(c322) <--> null
%		ks322:=  ;  kdeg322:=  ;
%	  % XIAP(c108) <--> null
%		ks108:=  ;  kdeg108:=  ;
%	  % PARP(c105) <--> null
%		ks105:=  ;  kdeg105:=  ;
%	  % Bid(c181) <--> null
%		ks181:=  ;  kdeg181:=  ;
%	  % Mcl1(c326) <--> null
%		ks326:=  ;  kdeg326:=  ;
%	  % Bax(c182) <--> null
%		ks182:=  ;  kdeg182:=  ;
%	  % Bcl2(c206) <--> null
%		ks206:=  ;  kdeg206:=  ;
%	  % M(c334) <--> null
%		ks334:=  ;  kdeg334:=  ;
%	  % CyCm(c188) <--> null
%		ks188:=  ;  kdeg188:=  ;
%	  % Smac(c339) <--> null
%		ks339:=  ;  kdeg339:=  ;
%	  % Apaf(c191) <--> null
%		ks191:=  ;  kdeg191:=  ;
%	  % C9(c193) <--> null
%		ks193:=  ;  kdeg193:=  ;
%	  % LR(c3) <--> null
%		ks3:= 0 ;  kdeg3:=  ;
%	  % R*(c317) <--> null
%		ks317:=  ;  kdeg317:=  ;
%	  % flip_R*(c318) <--> null
%		ks318:= 0  ;  kdeg318:=  ;
%	  % C8_R*(c319) <--> null
%		ks319:= 0 ;  kdeg319:=  ;
%	  % C8*(c56) <--> null
%		ks56:= 0 ;  kdeg563:=  ;
%	  % C8*_Bar(c321) <--> null
%		ks321:= 0 ;  kdeg321:=  ;
%	  % C8*_C3(c58) <--> null
%		ks58:= 0 ;  kdeg58:=  ;
%	  % C3*(c59) <--> null
%		ks59:= 0 ;  kdeg59:=  ;
%	  % C3*_C6(c323) <--> null
%		ks323:= 0 ;  kdeg323:=  ;
%	  % C6*(c324) <--> null
%		ks324:= 0 ;  kdeg324:=  ;
%	  % C6*_C8(c325) <--> null
%		ks325:= 0 ;  kdeg325:=  ;
%	  % XIAP_C3*(c210) <--> null
%		ks210:= 0 ;  kdeg210:=  ;
%	  % C3*_PARP(c106) <--> null
%		ks106:= 0 ;  kdeg106:=  ;
%	  % cPARP(c107) <--> cPARPdummy(c354)
%		ks107:= 0 ;  kdeg107:=  ;
%	  % C8*_Bid(c183) <--> null
%		ks183:= 0 ;  kdeg183:=  ;
%	  % tBid(c184) <--> null
%		ks184:= 0 ;  kdeg184:=  ;
%	  % Mcl1_tBid(c327) <--> null
%		ks327:= 0 ;  kdeg327:=  ;
%	  % Bax(c182) <--> null
%		ks182:= 0 ;  kdeg182:=  ;
%	  % tBid_Bax(c185) <--> null
%		ks185:= 0 ;  kdeg185:=  ;
%	  % Bax*(c328) <--> null
%		ks328:= 0 ;  kdeg328:=  ;
%	  % Bax*_m(c329) <--> null
%		ks329:= 0 ;  kdeg329:=  ;
%	  % Bax*_m_Bcl2(c330) <--> null
%		ks330:= 0 ;  kdeg330:=  ;
%	  % Bax2(c187) <--> null
%		ks187:= 0 ;  kdeg187:=  ;
%	  % Bax2_Bcl2(c331) <--> null
%		ks331:= 0 ;  kdeg331:=  ;
%	  % Bax4(c332) <--> null
%		ks332:= 0 ;  kdeg332:=  ;
%	  % Bax4_Bcl2(c333) <--> null
%		ks333:= 0 ;  kdeg333:=  ;
%	  % Bax4_M(c335) <--> null
%		ks335:= 0 ;  kdeg335:=  ;
%	  % M*(c336) <--> M
%		ks336:= 0 ;  kdeg336:=  ;
%	  % M*_Smac_m(c337) <--> null
%		ks337:= 0 ;  kdeg337:=  ;
%	  % CyCr(c338) <--> null
%		ks338:= 0 ;  kdeg338:=  ;
%	  % M*_Smac_m(c340) <--> null
%		ks340:= 0 ;  kdeg340:=  ;
%	  % Smac_r(c341) <--> null
%		ks341:= 0 ;  kdeg341:=  ;
%	  % CyC(c190) <--> null
%		ks190:= 0 ;  kdeg190:=  ;
%	  % Apaf_CyC(c192) <--> null
%		ks192:= 0 ;  kdeg192:=  ;
%	  % Apaf*(c342) <--> null
%		ks342:= 0 ;  kdeg342:=  ;
%	  % Apop(c343) <--> null
%		ks343:= 0 ;  kdeg343:=  ;
%	  % Apop_C3(c344) <--> null
%		ks344:= 0 ;  kdeg344:=  ;
%	  % Smac(c345) <--> SmacDummy(c355)
%		ks345:= 0 ;  kdeg345:=  ;
%	  % Apop_XIAP(c346) <--> null
%		ks346:= 0 ;  kdeg346:=  ;
%	  % Smac_XIAP(c347) <--> SmacDummy(c355)
%		ks347:= 0 ;  kdeg347:=  ;
%	  % C3*_Ub(c223) <--> null
%		ks223:= 0 ; kdeg223:=  ;
%	  % IETD(c348) <--> null
%		ks348:= 0 ; kdeg348:=  ;
%	  % C8:IETD(c349) <--> null
%		ks349:= 0 ; kdeg349:=  ;
%	  % DEVDR(c350) <--> null
%		ks350:= 0 ; kdeg350:=  ;
%	  % C3:DEVDR(c351) <--> null
%		ks351:= 0 ; kdeg351:=  ;
%		
