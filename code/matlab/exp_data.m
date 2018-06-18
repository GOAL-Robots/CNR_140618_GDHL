
function gadata = exp_data(np,idx,expname)
addpath('./expd')

%%
% BiPoo2001                         original
% CassenaerLaurent2007              extrapolated
% FroemkeDan2002                    extrapolated
% WittenbergWang20062               extrapolated
% WoodinGangulyPoo2003              original
% NishiyamaHongKato2000             original
% TzounopoulosKimTrussell20042      extrapolated
% ZhangTaoPoo1998                   extrapolated
% ZhouAckerWhite2005                original
% HaasNowotnyAbarbanel2006          original
%%

gadata.exp_name = expname; % AUTHORS
gadata.perm_idx = idx;
gadata.n_ex_par  = 3;
gadata.n_in_par  = np;
gadata.num_gen  = 100;
gadata.num_pop  = 600;
gadata.num_var  = gadata.n_ex_par + gadata.n_in_par;
gadata.lambda = 0.0000;

%gadata.crossover  = .01;
%gadata.penalty  = 10;
%gadata.scale_gap   = [ 0 30,  30, 100,   2,   2,   2,   2,   2,   2,   2,   2 ];

gadata.crossover  = .1;
gadata.penalty  = 1;
gadata.scale_gap   = [ 0 30,  30, 200,   2,   2,   2,   2,   2,   2,   2,   2 ];

gadata.scale_b     = [ 0   0,   0,   0,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1 ];


% permutations
p = perms(1:8);
p = unique(max(0,p.*(p>(8-gadata.n_in_par))-(8-gadata.n_in_par)),'rows');
p = unique(p>0,'rows');

p = p(gadata.perm_idx,:);
gadata.p = p;


c=1; dd=[]; d=[1 1 1 p];
for x = 1:length(d), 
  if d(x)>0,
    dd=[dd c];
    c=c+1;
  else
    dd=[dd 0];
  end;
end
gadata.pidx = dd;
switch gadata.exp_name
case 'BiPoo2001'
        load BiPoo2001;
        gadata.scal_f = 1.0;
        BiPoo2001fil = BiPoo2001;
        gadata.stdp_time = BiPoo2001fil(:,1);
        gadata.stdp_fun = BiPoo2001fil(:,2)./gadata.scal_f;
    case 'CassenaerLaurent2007'
        load CassenaerLaurent2007;
        gadata.scal_f = 1;
        gadata.stdp_time = CassenaerLaurent2007(:,1);
        gadata.stdp_fun = CassenaerLaurent2007(:,2)./gadata.scal_f;
    case 'FroemkeDan2002'
        load FroemkeDan2002;        
        %FroemkeDan2002 =stdp_data_bar(FroemkeDan2002,20);
        gadata.scal_f = 1;
        gadata.stdp_time = FroemkeDan2002(:,1);
        gadata.stdp_fun = FroemkeDan2002(:,2)./gadata.scal_f;
    case 'WittenbergWang20062'
        load WittenbergWang20062;
        WittenbergWang20062 = sortrows(WittenbergWang20062);
        %WittenbergWang20062 =stdp_data_bar(WittenbergWang20062,40);
        gadata.scal_f = 1;
        gadata.stdp_time = WittenbergWang20062(:,1);
        gadata.stdp_fun = (WittenbergWang20062(:,2)*50-30)./gadata.scal_f;
    case 'WoodinGangulyPoo2003'
        load WoodinGangulyPoo2003;
        gadata.scal_f = 1.0;
        gadata.stdp_time = WoodinGangulyPoo2003(:,1);
        gadata.stdp_fun = WoodinGangulyPoo2003(:,2)./gadata.scal_f;
    case 'NishiyamaHongKato2000'
        load NishiyamaHongKato2000;
        gadata.scal_f = 1.0;
        gadata.stdp_time = NishiyamaHongKato2000(:,1);
        gadata.stdp_fun = NishiyamaHongKato2000(:,2)./gadata.scal_f;
    case 'TzounopoulosKimTrussell20042'
        load TzounopoulosKimTrussell20042;
        gadata.scal_f = 1.0;
        gadata.stdp_time = TzounopoulosKimTrussell20042(:,1);
        gadata.stdp_fun = TzounopoulosKimTrussell20042(:,2)./gadata.scal_f;
    case 'ZhangTaoPoo1998'
        load ZhangTaoPoo1998;
        gadata.scal_factor = [10 10000 1000];
        gadata.scal_f = 1.0;
        gadata.stdp_time = ZhangTaoPoo1998(:,1);
        gadata.stdp_fun = ZhangTaoPoo1998(:,2)./gadata.scal_f;
    case 'ZhouAckerWhite2005'
        load ZhouAckerWhite2005;
        gadata.scal_f = 1.0;
        gadata.stdp_time = ZhouAckerWhite2005(:,1);
        gadata.stdp_fun = ZhouAckerWhite2005(:,2)./gadata.scal_f;
    case 'ZhouAckerWhite2005broad'
        load ZhouAckerWhite2005broad;
        gadata.scal_f = 1.0;
        gadata.stdp_time = ZhouAckerWhite2005broad(:,1);
        gadata.stdp_fun = ZhouAckerWhite2005broad(:,2)./gadata.scal_f;
    case 'HaasNowotnyAbarbanel2006'
        load HaasNowotnyAbarbanel2006;
        gadata.scal_f = 0.01;
        gadata.stdp_time = HaasNowotnyAbarbanel2006(:,1);
        gadata.stdp_fun = HaasNowotnyAbarbanel2006(:,2)./gadata.scal_f;
end
