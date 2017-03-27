k = 1;

for h = 1:length(epsilon_refined)
    
    j = 1;
    
    for i = 1:length(acceptable)
        
        if acceptable(1,i) == epsilon_refined(h)
            
            sig_scan(j) = acceptable(2,i);
            
            j = j + 1;
                        
        end
        
    end
    
    if j > 1
        
        sig_low(k) = min(sig_scan);
        sig_high(k) = max(sig_scan);
        eps_acc(k) = epsilon_refined(h);
       
        k = k+1;
        
    end
    
    clear sig_scan
    
end



% index = 1;
% 
% sig_low_index = 1;
% sig_high_index = 1;
% j = 1;
% 
% for i = 1:length(acceptable)
%    
%     eps_acceptable(j) = acceptable(i);
%     
%     
%     
%     
% end