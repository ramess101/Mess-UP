
s = 1;

for i = 1:length(sigma)-1
    
    for j = i+1:length(sigma)
        
        dif = abs(sigma(i) - sigma(j));
        
        if dif < 0.000002 
            
            dif2 = abs(epsilon(i) - epsilon(j));
            
            if dif2 < 0.00002
                
                duplicate(s) = i;
                               
                s = s + 1;
                
                duplicate(s) = j;
                
            end
            
        end
        
    end
    
end
        