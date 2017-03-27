clear

b_inherent = [0.206, 305.3, 0.000557 0.05557];
beta = 0.32;

T_data = linspace(175,300,50);

n = length(T_data);
p = 4;

Y_hat = @(b,T) b(1) + b(3) .* (b(2) - T) + b(4) .* (b(2) - T).^beta;

y_inherent = Y_hat(b_inherent,T_data);

TC_inherent = b_inherent(2);

s_inherent = 8.27*10^-4; % Constant uncertainty

valid = 0;
valid_TC = 0;

for a= 1:1
    
    b_guess = b_inherent;
    
    Y = normrnd(y_inherent,s_inherent);
        
    SSE_hat = @(b) sum((Y_hat(b,T_data) - Y).^2);
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    Y_fit = Y_hat(b_fit,T_data);
    
    SSE_fit = SSE_hat(b_fit);
    
    sigma2 = SSE_fit / (n - p);
    
    RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));
    
    % This is good for 20 temperatures between 175 and 275 with 8.27x10^-4
    % error
%     rhoc_range = linspace(0.7,1.3,30)*b_fit(1);
%     TC_range = linspace(0.9,1.1,30)*b_fit(2);
%     A_range = linspace(0.75,1.25,30)*b_fit(3);
%     B_range = linspace(0.7,1.3,30)*b_fit(4);
    
    % This is good for 50 temperatures between 175 and 300 with 8.27x10^-4
    % error
    rhoc_range = linspace(0.9,1.1,30)*b_fit(1);
    TC_range = linspace(0.99,1.01,30)*b_fit(2);
    A_range = linspace(0.85,1.15,30)*b_fit(3);
    B_range = linspace(0.85,1.15,30)*b_fit(4);

    k = 1;

    for g = 1:length(rhoc_range)
    
        for h=1:length(TC_range)
            
            for i = 1:length(A_range)
                
                for j = 1:length(B_range)
        
            b_range = [rhoc_range(g); TC_range(h); A_range(i); B_range(j)];
                    
            SSE_range = SSE_hat(b_range);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_hat(b_range,T_data);
                
                    k = k+1;
                    
                    if g == 1 || g == length(rhoc_range) || h == 1 || h == length(TC_range) || i == 1 || i == length(A_range) || j == 1 || j == length(B_range)
                       
                        pause(0)
                        
                    end
                                             
                end
                
                end
                
            end
            
        end
    
    end
    
    Y_hi = max(Y_ext);
    Y_lo = min(Y_ext);
    TC_hi = max(b_ext(:,2));
    TC_lo = min(b_ext(:,2));
        
    valid_all = 1;
           
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < y_inherent(r) || Y_lo(r) > y_inherent(r)
          
           valid_all = 0;
           
       end
       
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
    end
    
    if TC_hi >= TC_inherent && TC_lo <= TC_inherent
        
        valid_TC = valid_TC + 1;
        
    end
    
    figure
    hold
    subplot(3,2,1)
    plot(b_ext(:,1),b_ext(:,2))
    subplot(3,2,2)
    plot(b_ext(:,1),b_ext(:,3))
    subplot(3,2,3)
    plot(b_ext(:,1),b_ext(:,4))
    subplot(3,2,4)
    plot(b_ext(:,2),b_ext(:,3))
    subplot(3,2,5)
    plot(b_ext(:,2),b_ext(:,4))
    subplot(3,2,6)
    plot(b_ext(:,3),b_ext(:,4))
    
    max(b_ext(:,1))/b_fit(1)
    max(b_ext(:,2))/b_fit(2)
    max(b_ext(:,3))/b_fit(3)
    max(b_ext(:,4))/b_fit(4)

    clear b_ext Y_ext
    
end

valid/a
valid_TC/a

figure
hold
plot(T_data,y_inherent,'g')
scatter(T_data,Y)
plot(T_data,Y_fit)
plot(T_data,Y_hi,'--')
plot(T_data,Y_lo,'--')
hold


