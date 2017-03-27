eps_data_temp = zeros(1,10);

s = 1;

for i = 1:length(epsilon)
    
    for j = 1:length(eps_ext)
        
        if epsilon(i) == eps_ext(j) && sigma(i) == sig_ext(j)
            
            s = s + 20;
            
            break
            
        end
            
    end       
    
        eps_data_20 = eps_data(s:s+19,:);
            
        eps_data_temp = [eps_data_temp; eps_data_20];
        
        s = s + 20;
    
end

eps_data = eps_data_temp(2:end,:);