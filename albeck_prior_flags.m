function prior_flag=albeck_prior_flags()
%not a rate contant: 0
%forward, 1st order: 1
%reverse, 1st order: 2
%foward, 2nd order: 3
%reverse, 2nd order: 4
%catalytic: 5

prior_flag(1:31)=3;
prior_flag(13)=1;
prior_flag(22)=1;
prior_flag(26)=1;
prior_flag(31)=1;
prior_flag(32:62)=2;
prior_flag(63:78)=5;

end