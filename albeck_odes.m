function out = albeck_odes(t, conc, position)
 kf=position(1:31);
 kr=position(32:62);
 n=63;
for m=[1,3,5:10,12,19:21, 23, 25, 29, 30]
        kc(m)=position(n);
        n=n+1;
 end;
 
global  v kdeg ks
%list of ODEs for S. Gaudet's TRAIL model

% #============================================================================
% #    OpenBio model
% #    Model category: signal pathway
% #    Model title: caspase_cdp
% #    Model author/reference: Albeck et al. (2007)
% #    OpenBio model developer: T. Park (2007)
% #    Revision History
% #    v1		tsp		made OpenBio compliant
% #	 Edited by S. Gaudet 
% #	 v2		sg		Added synthesis and degratation terms with new parameter values
% #	 v3		sg		Changed species name to make coherent with TNF model
% #============================================================================

%this m-file is used by jacobian.m
    
    out(1,1) = -kf(1)*conc(1)*conc(2) +kr(1)*conc(3) +kf(31)*conc(4) -kdeg(1)*conc(1) +ks(1) ; % L(c2)  ks, kdeg(1) => 2  
    out(2,1) = -kf(1)*conc(2)*conc(1) +kr(1)*conc(3) +kf(31)*conc(4) -kdeg(2)*conc(2) +ks(2) ; % R(c1) ks, kdeg(2) => 1
    out(3,1) =  kf(1)*conc(1)*conc(2) -kr(1)*conc(3) -kc(1)*conc(3) -kdeg(3)*conc(3);% +ks(3); L:R(c3) ks, kdeg(3) => 3
    
    out(4,1) =  kc(1)*conc(3) -kf(2)*conc(4)*conc(5) +kr(2)*conc(6) ...
            -kf(3)*conc(4)*conc(7) +kr(3)*conc(8) +kc(3)*conc(8) -kf(31)*conc(4) -kdeg(4)*conc(4); % +ks(4); % R*(c317) ks, kdeg(4) => 317
   
    out(5,1) = -kf(2)*conc(4)*conc(5) +kr(2)*conc(6) -kdeg(5)*conc(5) +ks(5); % flip(c35) ks, kdeg(5) => 35
    out(6,1) =  kf(2)*conc(4)*conc(5) -kr(2)*conc(6) -kdeg(6)*conc(6); % +ks(6); % flip:R*(c318) ks, kdeg(6) => 318
    out(7,1) = -kf(3)*conc(4)*conc(7) +kr(3)*conc(8) ...
            -kf(7)*conc(7)*conc(17) +kr(7)*conc(18) -kdeg(7)*conc(7) +ks(7); % C8(c27) ks, kdeg(7) => 27
    out(8,1) =  kf(3)*conc(4)*conc(7) -kr(3)*conc(8) -kc(3)*conc(8) -kdeg(8)*conc(8); % +ks(8); % C8:R*(c319) ks, kdeg(8) => 319
    out(9,1) =  kc(3)*conc(8) -kf(4)*conc(9)*conc(10) +kr(4)*conc(11) ...
            -kf(5)*conc(9)*conc(12) +kr(5)*conc(13) +kc(5)*conc(13) +kc(7)*conc(18)... 
            -kf(10)*conc(9)*conc(24) +kr(10)*conc(25) +kc(10)*conc(25) ...
            -kf(29)*conc(9)*conc(59) +kr(29)*conc(60) +kc(29)*conc(60) -kdeg(9)*conc(9); % +ks(9); % C8*(c56) ks, kdeg(9) => 56
    out(10,1) = -kf(4)*conc(9)*conc(10) +kr(4)*conc(11) -kdeg(10)*conc(10) +ks(10); % Bar(c320) ks, kdeg(10) => 320
    out(11,1) =  kf(4)*conc(9)*conc(10) -kr(4)*conc(11) -kdeg(11)*conc(11); % +ks(11); % C8*:Bar(c321) ks, kdeg(11) => 321
    out(12,1) = -kf(5)*conc(9)*conc(12) +kr(5)*conc(13) .....
             -kf(25)*conc(12)*conc(53) +kr(25)*conc(54) -kdeg(12)*conc(12) +ks(12); % C3(c57) ks, kdeg(12) => 57
    out(13,1) =  kf(5)*conc(9)*conc(12) -kr(5)*conc(13) -kc(5)*conc(13) -kdeg(13)*conc(13); % +ks(13); % C8*:C3(c58) ks, kdeg(13) => 58
    out(14,1) =  kc(5)*conc(13) .....
             -kf(6)*conc(14)*conc(15) +kr(6)*conc(16) +kc(6)*conc(16) ...
             -kf(8)*conc(14)*conc(19) +kr(8)*conc(20) ...
             -kf(9)*conc(14)*conc(21) +kr(9)*conc(22) +kc(9)*conc(22) ...
             -kf(30)*conc(14)*conc(62) +kr(30)*conc(63) +kc(30)*conc(63)...
             +kc(25)*conc(54) -kdeg(14)*conc(14); % +ks(14) ; % C3*(c59)  ks, kdeg(14) => 59
    out(15,1) = -kf(6)*conc(14)*conc(15) +kr(6)*conc(16) -kdeg(15)*conc(15) +ks(15); % C6(c322)  ks, kdeg(15) => 322
    out(16,1) =  kf(6)*conc(14)*conc(15) -kr(6)*conc(16) -kc(6)*conc(16)-kdeg(16)*conc(16); % +ks(16) ; % C3*:C6(c323)  ks, kdeg(16) => 323
    out(17,1) =  kc(6)*conc(16) ...
             -kf(7)*conc(7)*conc(17) +kr(7)*conc(18) +kc(7)*conc(18)-kdeg(17)*conc(17); % +ks(17) ; % C6*(c324)  ks, kdeg(17) => 324
    out(18,1) =  kf(7)*conc(7)*conc(17) -kr(7)*conc(18) -kc(7)*conc(18) -kdeg(18)*conc(18); % +ks(18); % C6*:C8(c325)  ks, kdeg(18) => 325
    out(19,1) = -kf(8)*conc(14)*conc(19) +kr(8)*conc(20) +kc(8)*conc(20) ...
             -kf(27)*conc(19)*conc(53) +kr(27)*conc(56) ...
             -kf(28)*conc(19)*conc(55) +kr(28)*conc(57) -kdeg(19)*conc(19) +ks(19); % XIAP(c108)  ks, kdeg(19) => 108
    out(20,1) =  kf(8)*conc(14)*conc(19) -kr(8)*conc(20) -kc(8)*conc(20) -kdeg(20)*conc(20); % +ks(20); % xiap:C3*(c210)  ks, kdeg(20) => 210
    out(21,1) = -kf(9)*conc(14)*conc(21) +kr(9)*conc(22) -kdeg(21)*conc(21) +ks(21); % PARP(c105)  ks, kdeg(21) => 105
    out(22,1) =  kf(9)*conc(14)*conc(21) -kr(9)*conc(22) -kc(9)*conc(22) -kdeg(22)*conc(22); % +ks(22); % C3*:PARP(c106)  ks, kdeg(22) => 106
    out(23,1) = kc(9)*conc(22) -kdeg(23)*conc(23); % +ks(23); % cPARP(c107)  ks, kdeg(23) => 107
    out(24,1) = -kf(10)*conc(9)*conc(24) +kr(10)*conc(25) -kdeg(24)*conc(24) +ks(24); % Bid(c181)  ks, kdeg(24) => 181
    out(25,1) =  kf(10)*conc(9)*conc(24) -kr(10)*conc(25) -kc(10)*conc(25) -kdeg(25)*conc(25); % +ks(25); % C8*:Bid(c183)  ks, kdeg(25) => 183
    out(26,1) =  kc(10)*conc(25) -kf(11)*conc(26)*conc(27) +kr(11)*conc(28) .....
             -kf(12)*conc(26)*conc(29) +kr(12)*conc(30) + kc(12)*conc(30)-kdeg(26)*conc(26); % +ks(26); % tBid(c184)  ks, kdeg(26) => 184
    out(27,1) = -kf(11)*conc(26)*conc(27) +kr(11)*conc(28) -kdeg(27)*conc(27) +ks(27); % Mcl1(c326)  ks, kdeg(27) => 326
    out(28,1) = +kf(11)*conc(26)*conc(27) -kr(11)*conc(28)-kdeg(28)*conc(28); % +ks(28); % tBid:Mcl1(c327)  ks, kdeg(28) => 327
    out(29,1) = -kf(12)*conc(26)*conc(29) +kr(12)*conc(30) -kdeg(29)*conc(29) +ks(29); % Bax(c182)  ks, kdeg(29) => 182
    out(30,1) =  kf(12)*conc(26)*conc(29) -kr(12)*conc(30) - kc(12)*conc(30) -kdeg(30)*conc(30); % +ks(30); % tBid:Bax(c185)  ks, kdeg(30) => 185
    out(31,1) =  kc(12)*conc(30) ...
             -kf(13)*conc(31) + kr(13)*conc(32) -kdeg(31)*conc(31); % +ks(31); % Bax*(c328)  ks, kdeg(31) => 328
    out(32,1) =  kf(13)*conc(31) - kr(13)*conc(32) -1/v*kf(14)*conc(32)*conc(33) +kr(14)*conc(34) ...
             -1/v*2*kf(15)*conc(32)*conc(32) +2*kr(15)*conc(35) -kdeg(32)*conc(32); % +ks(32); % Bax*_m(c329)  ks, kdeg(32) => 329
    out(33,1) = -1/v*kf(14)*conc(32)*conc(33) +kr(14)*conc(34) ...
             -1/v*kf(16)*conc(33)*conc(35) +kr(16)*conc(36) ...
             -1/v*kf(18)*conc(33)*conc(37) +kr(18)*conc(38) -kdeg(33)*conc(33) +ks(33); % Bcl2(c206)  ks, kdeg(33) => 206
    out(34,1) =  1/v*kf(14)*conc(32)*conc(33) -kr(14)*conc(34) -kdeg(34)*conc(34); % +ks(34); % Bax*_m:Bcl2(c330)  ks, kdeg(34) => 330
    out(35,1) =  1/v*kf(15)*conc(32)*conc(32) -kr(15)*conc(35)...
             -1/v*kf(16)*conc(33)*conc(35) +kr(16)*conc(36) ...
             -2/v*kf(17)*conc(35)*conc(35) +2*kr(17)*conc(37) -kdeg(35)*conc(35); % +ks(35); % Bax2(c187)  ks, kdeg(35) => 187
    out(36,1) =  1/v*kf(16)*conc(33)*conc(35) -kr(16)*conc(36) -kdeg(36)*conc(36); % +ks(36); % Bax2:Bcl2(c331)  ks, kdeg(36) => 331
    out(37,1) = 1/v*kf(17)*conc(35)*conc(35) -kr(17)*conc(37)...
            -1/v*kf(18)*conc(33)*conc(37) +kr(18)*conc(38) ...
            -1/v*kf(19)*conc(39)*conc(37) +kr(19)*conc(40) -kdeg(37)*conc(37); % +ks(37); % Bax4(c332)  ks, kdeg(37) => 332
    out(38,1) = 1/v*kf(18)*conc(33)*conc(37) -kr(18)*conc(38) -kdeg(38)*conc(38); % +ks(38); % Bax4:Bcl2(c333)  ks, kdeg(38) => 333
    out(39,1) = -1/v*kf(19)*conc(39)*conc(37) +kr(19)*conc(40) -kdeg(39)*conc(39)  +kdeg(41)*conc(41) +ks(39); %   M(c334)  ks, kdeg(39) => 334
    out(40,1) =  1/v*kf(19)*conc(39)*conc(37) -kr(19)*conc(40) -kc(19)*conc(40) -kdeg(40)*conc(40); % +ks(40); % Bax4:M(c335)  ks, kdeg(40) => 335
    out(41,1) =  kc(19)*conc(40) ...
             -1/v*kf(20)*conc(41)*conc(42) +kr(20)*conc(43) +kc(20)*conc(43) ...
             -1/v*kf(21)*conc(41)*conc(45) +kr(21)*conc(46) +kc(21)*conc(46) -kdeg(41)*conc(41); % +ks(41); % M*(c336)  ks, kdeg(41) => 336
    out(42,1) = -1/v*kf(20)*conc(41)*conc(42) +kr(20)*conc(43) -kdeg(42)*conc(42) +ks(42); % CyCm(c188)  ks, kdeg(42) => 188
    out(43,1) =  1/v*kf(20)*conc(41)*conc(42) -kr(20)*conc(43) -kc(20)*conc(43) -kdeg(43)*conc(43); % +ks(43); % M*:CyCm(c337)  ks, kdeg(43) => 337
    out(44,1) =  kc(20)*conc(43) ...
             -kf(22)*conc(44) +kr(22)*conc(48) -kdeg(44)*conc(44); % +ks(44); % CyCr(c338)  ks, kdeg(44) => 338
    out(45,1) = -1/v*kf(21)*conc(41)*conc(45) +kr(21)*conc(46) -kdeg(45)*conc(45) +ks(45); % Smac_m(c339)  ks, kdeg(45) => 339
    out(46,1) =  1/v*kf(21)*conc(41)*conc(45) -kr(21)*conc(46) -kc(21)*conc(46) -kdeg(46)*conc(46); % +ks(46); % M*:Smac_m(c340)  ks, kdeg(46) => 340
    out(47,1) =  kc(21)*conc(46) ...
             -kf(26)*conc(47) +kr(26)*conc(55) -kdeg(47)*conc(47); % +ks(47); % Smac_r(c341)  ks, kdeg(47) => 341
    out(48,1) =  kf(22)*conc(44) -kr(22)*conc(48) ...
             -kf(23)*conc(48)*conc(49) +kr(23)*conc(50) +kc(23)*conc(50) -kdeg(48)*conc(48); % +ks(48); % CyC(c190)  ks, kdeg(48) => 190
    out(49,1) = -kf(23)*conc(48)*conc(49) +kr(23)*conc(50) -kdeg(49)*conc(49) +ks(49); % Apaf(c191)  ks, kdeg(49) => 191
    out(50,1) =  kf(23)*conc(48)*conc(49) -kr(23)*conc(50) -kc(23)*conc(50) -kdeg(50)*conc(50); % +ks(50); % Apaf:CyC(c192)  ks, kdeg(50) => 192
    out(51,1) = kc(23)*conc(50)...
            -kf(24)*conc(51)*conc(52) +kr(24)*conc(53) -kdeg(51)*conc(51); % +ks(51); % Apaf*(c342)  ks, kdeg(51) => 342
    out(52,1) = -kf(24)*conc(51)*conc(52) +kr(24)*conc(53) -kdeg(52)*conc(52) +ks(52); % C9(c193)  ks, kdeg(52) => 193
    out(53,1) =  kf(24)*conc(51)*conc(52) -kr(24)*conc(53)...
             -kf(25)*conc(12)*conc(53) +kr(25)*conc(54) +kc(25)*conc(54)...
             -kf(27)*conc(19)*conc(53) +kr(27)*conc(56) -kdeg(53)*conc(53); % +ks(53); % Apop(c343)  ks, kdeg(53) => 343
    out(54,1) =  kf(25)*conc(12)*conc(53) -kr(25)*conc(54) -kc(25)*conc(54) -kdeg(54)*conc(54); % +ks(54); % Apop:C3(c344)  ks, kdeg(54) => 344
    out(55,1) =  kf(26)*conc(47) -kr(26)*conc(55)...
             -kf(28)*conc(19)*conc(55) +kr(28)*conc(57) -kdeg(55)*conc(55); % +ks(55); % Smac(c345)  ks, kdeg(55) => 345
    out(56,1) =  kf(27)*conc(19)*conc(53) -kr(27)*conc(56) -kdeg(56)*conc(56); % +ks(56); % Apop:XIAP(c346)  ks, kdeg(56) => 346
    out(57,1) =  kf(28)*conc(19)*conc(55) -kr(28)*conc(57) -kdeg(57)*conc(57); % +ks(57); % Smac:XIAP(c347)  ks, kdeg(57) => 347
    out(58,1) =  kc(8)*conc(20) -kdeg(58)*conc(58); % +ks(58); % C3*_Ub(c223)  ks, kdeg(58) => 223
	out(59,1) = -kf(29)*conc(59)*conc(9) +kr(29)*conc(60) -kdeg(59)*conc(59) +ks(59); % IETD(c348)  ks, kdeg(59) => 348
	out(60,1) =  kf(29)*conc(59)*conc(9) -kr(29)*conc(60) -kc(29)*conc(60) -kdeg(60)*conc(60); % +ks(60); % IETD:C8*(c349)  ks, kdeg(60) => 349
	out(61,1) =  kc(29)*conc(60) -kdeg(61)*conc(61); % +ks(61);  % IETD*(c350)  ks, kdeg(61) => 350
	out(62,1) = -kf(30)*conc(62)*conc(14) +kr(30)*conc(63) -kdeg(62)*conc(62) +ks(62); % DEDVR(c351)  ks, kdeg(62) => 351
	out(63,1) =  kf(30)*conc(62)*conc(14) -kr(30)*conc(63) -kc(30)*conc(63) -kdeg(63)*conc(63); % +ks(63); % DEDVR:C3*(c352)  ks, kdeg(63)  => 352
	out(64,1) =  kc(30)*conc(63) -kdeg(64)*conc(64); % +ks(64);  % DEDVR*(c353)  ks, kdeg(63) => 353
	out(65,1) =  kdeg(23)*conc(23);  % accumulation of a cPARPdummy species, to account for all created, c354
	out(66,1) =  kdeg(47)*conc(47) +kdeg(55)*conc(55) +kdeg(57)*conc(57) ;  % accumulation of a SMACdummy species, to account for all created, c355
	out(67,1) =  kc(9)*conc(22) -kdeg(23)*conc(23) +kdeg(23)*conc(23);  %c356
	out(68,1) =  kc(21)*conc(46) -kf(26)*conc(47) +kr(26)*conc(55) -kdeg(47)*conc(47)...
			+kf(26)*conc(47) -kr(26)*conc(55) -kf(28)*conc(19)*conc(55) +kr(28)*conc(57) -kdeg(55)*conc(55)...
			+kf(28)*conc(19)*conc(55) -kr(28)*conc(57) -kdeg(57)*conc(57) ; %c357
	out(69) = kc(21)*conc(46) -kf(26)*conc(47) +kr(26)*conc(55) -kdeg(47)*conc(47)...
			+kf(26)*conc(47) -kr(26)*conc(55) -kf(28)*conc(19)*conc(55) +kr(28)*conc(57) -kdeg(55)*conc(55)...
			+kf(28)*conc(19)*conc(55) -kr(28)*conc(57) -kdeg(57)*conc(57)...
			+kdeg(47)*conc(47) +kdeg(55)*conc(55) +kdeg(57)*conc(57) ;  %c358
        
end



     

	