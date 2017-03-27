clear

beta_range = linspace(0.3,0.33,100);

for bb = 1:length(beta_range)
    
    beta = beta_range(bb);

    Parameter_space_analysis_ethane_alternative

    for h = 1:length(epsilon)
    
        for i =1:length(sigma)
            
            if bb == 1
               
                SSE_opt(h,i) = SSE(h,i);
                
            elseif SSE(h,i) < SSE_opt(h,i)
            
                SSE_opt(h,i) = SSE(h,i);
                beta_opt(h,i) = beta;
                TC_opt(h,i) = TC(h,i);

            end
            
        end
    
    end

end