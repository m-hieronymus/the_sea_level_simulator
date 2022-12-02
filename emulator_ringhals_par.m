clear all; close all
tic
load distributions
nr_res=500;

get_liklihoodCI=0; % calculate distribution for GEV parameters
ext_par_uncert=0; %GEV parameters are not known if 1, but their distribution is;parpool([1 100]) 
nr_iter=1000; %iterations times nr_par
nr_par=10000; %number of parfor iterations
p19=0.05;
p26=0.155;
p26l=0.01;
p45=0.5;
p7=0.22;
p85=0.064;
p85l=1-p19-p26-p45-p7-p85-p26l;
if p85l<0
    disp('p85 cannot be < 0')
    return
end
T=[2020 2030 2040 2050 2060 2070 2080 2090 2100 2110 2120 2130 2140 2150];
t=2020:2150;
load amax_Ringhals
[m_paramEsts,m_parmci] = gevfit(amaxs);
    obs_kMLE = m_paramEsts(1);        % Shape parameter
    obs_sigmaMLE = m_paramEsts(2);    % Scale parameter
    obs_muMLE = m_paramEsts(3);       % Location parameter
x=11:10:140;
pX=length(x);
long_n=150;

ssh_ext=linspace(0.3,3.7,nr_res);
ssh_mean=linspace(0,12,nr_res);
ssh_joint=linspace(0,14,nr_res);


if get_liklihoodCI==1
    probs=0.01:0.01:0.99;
    parms_l=nan(99,3);
    parms_l(:,1)=obs_kMLE;
    parms_l(:,2)=obs_sigmaMLE;
    parms_l(:,3)=obs_muMLE;
    parms_h=parms_l;
    kk=linspace(m_parmci(1,1)*0.9,m_parmci(2,1)*1.1,long_n);
    sigsig=linspace(m_parmci(1,2)*0.9,m_parmci(2,2)*1.1,long_n);
    mumu=linspace(m_parmci(1,3)*0.9,m_parmci(2,3)*1.1,long_n);
    %Chi squared confidence bounds
   for cl=1:99 
       lowest_ret=1000;
       highest_ret=0;
    % params=m_paramEsts;
    
    o_nllCritVal = gevlike([obs_kMLE,obs_sigmaMLE,obs_muMLE],amaxs) + .5*chi2inv(probs(cl),1); %95%
    for k=1:long_n
        for sig=1:long_n
            for mu=1:long_n
                if gevlike([kk(k),sigsig(sig),mumu(mu)],amaxs)<o_nllCritVal && gevinv(1-1/1000,kk(k),sigsig(sig),mumu(mu))>highest_ret
                    parms_h(cl,:)=[kk(k),sigsig(sig),mumu(mu)];
                    highest_ret=gevinv(1-1/1000,kk(k),sigsig(sig),mumu(mu));
                end
                if gevlike([kk(k),sigsig(sig),mumu(mu)],amaxs)<o_nllCritVal && gevinv(1-1/1000,kk(k),sigsig(sig),mumu(mu))<lowest_ret
                    parms_l(cl,:)=[kk(k),sigsig(sig),mumu(mu)];
                    lowest_ret=gevinv(1-1/1000,kk(k),sigsig(sig),mumu(mu));

                end
            end
        end
    end
   end
   parms=[m_paramEsts;parms_l;parms_h];
   save ext_parm_set parms
else
    load ext_parm_set
end

ssh_reached_joint=zeros(nr_res,pX);     %Reduction Variables
ssh_reached_mean_only=zeros(nr_res,pX);
ssh_reached_mean=zeros(nr_res,nr_res,pX);
ssh_reached_ext_only=zeros(nr_res,pX);
ssh_reached_ext=zeros(nr_res,nr_res,pX);

scen_ssh=zeros(nr_res,pX);
mean_quant_ssh=zeros(nr_res,pX);
ext_quant_ssh=zeros(nr_res,pX);

mean_of_tot=nan(nr_res,pX);
ext_of_tot=nan(nr_res,pX);


for a=1:nr_iter
    ssh_reached_tmp=nan(nr_par,5,pX);
    probs_tmp=nan(nr_par,3);
    parfor aa=1:nr_par
       
    p=rand(3,1);
    probs_tmp(aa,:)=p;
    if p(1)>0 && p(1)<p19
        mean_sea_coarse=[0,skn19_1.InverseCDF(p(2)),skn19_2.InverseCDF(p(2))...
            ,skn19_3.InverseCDF(p(2)),skn19_4.InverseCDF(p(2)),skn19_5.InverseCDF(p(2))...
            , skn19_6.InverseCDF(p(2)), skn19_7.InverseCDF(p(2)), skn19_8.InverseCDF(p(2)), skn19_9.InverseCDF(p(2)),  skn19_10.InverseCDF(p(2))...
            , skn19_11.InverseCDF(p(2)), skn19_12.InverseCDF(p(2)),skn19_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
    elseif p(1)>p19 && p(1)<p19+p26
        mean_sea_coarse=[0,skn26_1.InverseCDF(p(2)),skn26_2.InverseCDF(p(2))...
            ,skn26_3.InverseCDF(p(2)),skn26_4.InverseCDF(p(2)),skn26_5.InverseCDF(p(2))...
            , skn26_6.InverseCDF(p(2)), skn26_7.InverseCDF(p(2)), skn26_8.InverseCDF(p(2)), skn26_9.InverseCDF(p(2)),  skn26_10.InverseCDF(p(2))...
            , skn26_11.InverseCDF(p(2)), skn26_12.InverseCDF(p(2)),skn26_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    elseif p(1)>p19+p26 && p(1)<p19+p26+p26l
        mean_sea_coarse=[0,skn26low_1.InverseCDF(p(2)),skn26low_2.InverseCDF(p(2))...
            ,skn26low_3.InverseCDF(p(2)),skn26low_4.InverseCDF(p(2)),skn26low_5.InverseCDF(p(2))...
            , skn26low_6.InverseCDF(p(2)), skn26low_7.InverseCDF(p(2)), skn26low_8.InverseCDF(p(2)), skn26low_9.InverseCDF(p(2)),  skn26low_10.InverseCDF(p(2))...
            , skn26low_11.InverseCDF(p(2)), skn26low_12.InverseCDF(p(2)),skn26low_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    elseif p(1)>p19+p26+p26l && p(1)<p19+p26+p26l+p45
        mean_sea_coarse=[0,skn45_1.InverseCDF(p(2)),skn45_2.InverseCDF(p(2))...
            ,skn45_3.InverseCDF(p(2)),skn45_4.InverseCDF(p(2)),skn45_5.InverseCDF(p(2))...
            , skn45_6.InverseCDF(p(2)), skn45_7.InverseCDF(p(2)), skn45_8.InverseCDF(p(2)), skn45_9.InverseCDF(p(2)),  skn45_10.InverseCDF(p(2))...
            , skn45_11.InverseCDF(p(2)), skn45_12.InverseCDF(p(2)),skn45_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    elseif p(1)>p19+p26+p26l+p45 && p(1)<p19+p26+p26l+p45+p7
        mean_sea_coarse=[0,skn7_1.InverseCDF(p(2)),skn7_2.InverseCDF(p(2))...
            ,skn7_3.InverseCDF(p(2)),skn7_4.InverseCDF(p(2)),skn7_5.InverseCDF(p(2))...
            , skn7_6.InverseCDF(p(2)), skn7_7.InverseCDF(p(2)), skn7_8.InverseCDF(p(2)), skn7_9.InverseCDF(p(2)),  skn7_10.InverseCDF(p(2))...
            , skn7_11.InverseCDF(p(2)), skn7_12.InverseCDF(p(2)),skn7_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    elseif p(1)>p19+p26+p26l+p45+p7 && p(1)<p19+p26+p26l+p45+p7+p85
        mean_sea_coarse=[0,skn85_1.InverseCDF(p(2)),skn85_2.InverseCDF(p(2))...
            ,skn85_3.InverseCDF(p(2)),skn85_4.InverseCDF(p(2)),skn85_5.InverseCDF(p(2))...
            , skn85_6.InverseCDF(p(2)), skn85_7.InverseCDF(p(2)), skn85_8.InverseCDF(p(2)), skn85_9.InverseCDF(p(2)),  skn85_10.InverseCDF(p(2))...
            , skn85_11.InverseCDF(p(2)), skn85_12.InverseCDF(p(2)),skn85_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    elseif p(1)>p19+p26+p26l+p45+p7+p85 && p(1)<1
        mean_sea_coarse=[0,skn85low_1.InverseCDF(p(2)),skn85low_2.InverseCDF(p(2))...
            ,skn85low_3.InverseCDF(p(2)),skn85low_4.InverseCDF(p(2)),skn85low_5.InverseCDF(p(2))...
            , skn85low_6.InverseCDF(p(2)), skn85low_7.InverseCDF(p(2)), skn85low_8.InverseCDF(p(2)), skn85low_9.InverseCDF(p(2)),  skn85low_10.InverseCDF(p(2))...
            , skn85low_11.InverseCDF(p(2)), skn85low_12.InverseCDF(p(2)),skn85low_13.InverseCDF(p(2))];
        mean_sea_coarse(1)=mean_sea_coarse(2)-(mean_sea_coarse(3)-mean_sea_coarse(2));
        mean_sea_coarse=mean_sea_coarse-mean_sea_coarse(1);
        
    end
    mean_sea=interp1(T,mean_sea_coarse,t);
    
    if ext_par_uncert==0
        ext_sea=gevrnd(obs_kMLE,obs_sigmaMLE,obs_muMLE,1,length(t));
    else
        np=round(p(3)*198)+1;
        ext_sea=gevrnd(parms(np,1),parms(np,2),parms(np,3),1,length(t));
    end
   
    
    for b=1:pX
        
        
        [max_j,ind_jm]=max(ext_sea(1:x(b))+mean_sea(1:x(b)));
        [dummy,ind_j]=min(abs(max_j-ssh_joint));
        [dummy,ind_m]=min(abs(mean_sea(ind_jm)-ssh_mean));
        [dummy,ind_e]=min(abs(ext_sea(ind_jm)-ssh_ext));
        
        [max_e,dummy]=max(ext_sea(1:x(b)));
        [dummy,ind_e2]=min(abs(max_e-ssh_ext));
        
        [max_m,dummy]=max(mean_sea(1:x(b)));
        [dummy,ind_m2]=min(abs(max_m-ssh_mean));
        
        ssh_reached_tmp(aa,:,b)=[ind_j,ind_m,ind_e,ind_e2,ind_m2]
    end
    
    end
    for aaa=1:nr_par
        
        for bb=1:pX
           ssh_reached_joint(ssh_reached_tmp(aaa,1,bb),bb)=ssh_reached_joint(ssh_reached_tmp(aaa,1,bb),bb)+1; 
           ssh_reached_mean(ssh_reached_tmp(aaa,1,bb),ssh_reached_tmp(aaa,2,bb),bb)=ssh_reached_mean(ssh_reached_tmp(aaa,1,bb),ssh_reached_tmp(aaa,2,bb),bb)+1;
           ssh_reached_ext(ssh_reached_tmp(aaa,1,bb),ssh_reached_tmp(aaa,3,bb),bb)=ssh_reached_ext(ssh_reached_tmp(aaa,1,bb),ssh_reached_tmp(aaa,3,bb),bb)+1;
           ssh_reached_ext_only(ssh_reached_tmp(aaa,4,bb),bb)=ssh_reached_ext_only(ssh_reached_tmp(aaa,4,bb),bb)+1; 
           ssh_reached_mean_only(ssh_reached_tmp(aaa,5,bb),bb)=ssh_reached_mean_only(ssh_reached_tmp(aaa,5,bb),bb)+1; 

           scen_ssh(ssh_reached_tmp(aaa,1,bb),bb)=scen_ssh(ssh_reached_tmp(aaa,1,bb),bb)+probs_tmp(aaa,1);
           mean_quant_ssh(ssh_reached_tmp(aaa,1,bb),bb)=mean_quant_ssh(ssh_reached_tmp(aaa,1,bb),bb)+probs_tmp(aaa,2);
           ext_quant_ssh(ssh_reached_tmp(aaa,1,bb),bb)=ext_quant_ssh(ssh_reached_tmp(aaa,1,bb),bb)+probs_tmp(aaa,3);       
        end
    end
    
end
scen_ssh=scen_ssh./ssh_reached_joint;
mean_quant_ssh=mean_quant_ssh./ssh_reached_joint;
ext_quant_ssh=ext_quant_ssh./ssh_reached_joint;

cdf_ext_only_2150=1-cumsum(ssh_reached_ext_only(:,13))/(nr_iter*nr_par);
cdf_mean_only_2150=1-cumsum(ssh_reached_mean_only(:,13))/(nr_iter*nr_par);
cdf_joint_2150=1-cumsum(ssh_reached_joint(:,13))/(nr_iter*nr_par);

cdf_ext_only_2150(cdf_ext_only_2150==0)=nan;
cdf_mean_only_2150(cdf_mean_only_2150==0)=nan;
cdf_joint_2150(cdf_joint_2150==0)=nan;

cdf_ext_only_2050=1-cumsum(ssh_reached_ext_only(:,3))/(nr_iter*nr_par);
cdf_mean_only_2050=1-cumsum(ssh_reached_mean_only(:,3))/(nr_iter*nr_par);
cdf_joint_2050=1-cumsum(ssh_reached_joint(:,3))/(nr_iter*nr_par);

cdf_ext_only_2050(cdf_ext_only_2050==0)=nan;
cdf_mean_only_2050(cdf_mean_only_2050==0)=nan;
cdf_joint_2050(cdf_joint_2050==0)=nan;

cdf_ext_only_2100=1-cumsum(ssh_reached_ext_only(:,8))/(nr_iter*nr_par);
cdf_mean_only_2100=1-cumsum(ssh_reached_mean_only(:,8))/(nr_iter*nr_par);
cdf_joint_2100=1-cumsum(ssh_reached_joint(:,8))/(nr_iter*nr_par);

cdf_ext_only_2100(cdf_ext_only_2100==0)=nan;
cdf_mean_only_2100(cdf_mean_only_2100==0)=nan;
cdf_joint_2100(cdf_joint_2100==0)=nan;


ssh_reached_mean(ssh_reached_mean==0)=nan;
ssh_reached_ext(ssh_reached_ext==0)=nan;

for a=1:nr_res
    for b=1:pX
        mean_of_tot(a,b)=nansum(ssh_reached_mean(a,:,b).*ssh_mean)/nansum(ssh_reached_mean(a,:,b));
        ext_of_tot(a,b)=nansum(ssh_reached_ext(a,:,b).*ssh_ext)/nansum(ssh_reached_ext(a,:,b));
        
    end
end
toc




save example
