clear

T_Ar = [50.00 100.00 150.00 200.00 250.00 273.15 293.15 300.00 313.15 333.15 353.15 373.15 473.15 523.15 573.15 623.15 673.15 723.15 773.15 873.15 973.15 1073.15 1173.15 1273.15]; 
B_Ar = [-0.71346 -0.1908 -0.08575 -0.04848 -0.02877 -0.02241 -0.01782 -0.0164 -0.01387 -0.01045 -0.00744 -0.00479 0.0048 0.00805 0.01063 0.0127 0.01436 0.01572 0.01682 0.01847 0.01959 0.02036 0.02093 0.02137]; % L/mol 

% These are the values that are commonly reported in literature
sig_Ar = 3.405;
eps_Ar = 119;

T = linspace(50,500,30);

n = length(T);
p = 2;

B_inherent = second_virial(T,eps_Ar,sig_Ar);

s_inherent = 0.02 * mean(abs(B_inherent));

valid = 0;

for a= 1:1
    
    b_guess = [eps_Ar,sig_Ar];
    
    Y = normrnd(B_inherent,s_inherent);
    
%     If you want to see what the DIPPR uncertainty looks like
    T = T_Ar;
    Y = B_Ar; 
    
    Y_hat = @(b) second_virial(T,b(1),b(2));
    
%     SSE_hat = @(b) sum((Y_hat(b) - Y).^2);
    SSE_hat = @(b) sum(((Y_hat(b) - Y)./Y).^2); % This approach gives a result closer to what is reported
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    Y_fit = second_virial(T,b_fit(1),b_fit(2));
    
    SSE_fit = SSE_hat(b_fit);
    
    sigma2 = SSE_fit / (n - p);
    
    RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));
    
    % This is good for 30 Temperature values between 50 and 500 with 2%
    % error on the mean value
%     eps_range = linspace(0.975,1.025,100)*b_fit(1);
%     sig_range = linspace(0.975,1.025,100)*b_fit(2);
    
%     This is good for the actual DIPPR data from 50-373 without weighting
%     by relative error
%     eps_range = linspace(0.93,1.07,100)*b_fit(1);
%     sig_range = linspace(0.95,1.05,100)*b_fit(2);
    
    %     This is good for all of the DIPPR data with weighting by relative error
    eps_range = linspace(0.994,1.006,100)*b_fit(1);
    sig_range = linspace(0.993,1.007,100)*b_fit(2);


    k = 1;

    for i = 1:length(eps_range)
    
        for j=1:length(sig_range)
        
            b_range = [eps_range(i); sig_range(j)];
            
            Y_range = second_virial(T,b_range(1),b_range(2));
        
            SSE_range = SSE_hat(b_range);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_range;
                
                    k = k+1;
                    
                    if i == 1 || i == length(eps_range) || j == 1 || j == length(sig_range)
                       
                        pause(0)
                        
                    end
                                             
                end
            
        end
    
    end
    
    Y_hi = max(Y_ext);
    Y_lo = min(Y_ext);
        
    valid_all = 1;
        
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < B_inherent(r) || Y_lo(r) > B_inherent(r)
          
           valid_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
    end
    
    figure
    hold
    scatter(b_ext(:,2),b_ext(:,1))
    scatter(sig_Ar,eps_Ar,'g')
    hold
    
    max(b_ext(:,1))/b_fit(1)
    min(b_ext(:,1))/b_fit(1)
    max(b_ext(:,2))/b_fit(2)
    min(b_ext(:,2))/b_fit(2)
    
    clear Y_ext b_ext
    
end

valid/a

% file = fopen('B_95.txt','w');
% fwrite(file,valid/a)
% fclose(file);

figure
hold
plot(T,B_inherent,'g')
scatter(T,Y)
plot(T,Y_fit)
plot(T,Y_hi,'--')
plot(T,Y_lo,'--')
scatter(T_Ar,B_Ar)
hold


