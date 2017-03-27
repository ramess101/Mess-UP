clear rhoL_grapher eps_grapher sig_grapher

k = 1;

Temp_i = 4;

for i = 1:length(epsilon)
    
    for j = 1:length(sigma)
        
        rhoL_grapher(k,1) = liq_ave_ave_contour(i,j,Temp_i);
        
%         RMS_grapher(k,1) = RMS_opt(i,j);
        
        eps_grapher(k,1) = epsilon(i);
        
        sig_grapher(k,1) = sigma(j);
        
        k = k+1;
        
    end
    
end

k = 1;

for i = 1:length(epsilon_refined)
    
    for j = 1:length(sigma_refined)
        
        rhoL_grapher(k,2) = liq_hat_contour(i,j,Temp_i);
        
        RMS_grapher(k,1) = RMS_opt(i,j);
        
        eps_grapher(k,2) = epsilon_refined(i);
        
        sig_grapher(k,2) = sigma_refined(j);
        
        k = k+1;
        
    end
    
end

