% Need to run the MC file first to get all the vap_ave and liq_ave

clearvars -except Temp vap_ave liq_ave sigma

MW = 675.2878;

T = Temp;

PC = zeros(length(sigma),1);
ZC = PC;
sigma_PC = PC;
sigma_ZC = ZC;

for m = 61:length(sigma)
    
    pv = vap_ave(:,m)';
    pL = liq_ave(:,m)';
          
%     [TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);
%     
%     PC = PC_from_Rackett_D(TC_fit,pc_fit,MW/pc_fit,A_fit,b_fit,T(1)/TC_fit,2/7);
%        
%     ZC = PC * MW / (TC_fit * pc_fit * 8.314472);

    cd ../GEMC_Simulations_Analysis/
    
    [PC_low,PC_high,ZC_low,ZC_high] = rigorous_statistics_PC(T,pv,pL,MW);
    
    PC(m) = (PC_high + PC_low)/2;
    
    ZC(m) = (ZC_high + ZC_low)/2;
    
    sigma_PC(m) = (PC_high - PC_low) / 2 / 1.96;
    
    sigma_ZC(m) = (ZC_high - ZC_low) / 2 / 1.96;
    
    m
            
end