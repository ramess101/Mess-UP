clear

b_inherent = [2, 200];

x_data = linspace(50,500,30);

n = length(x_data);
p = 2;

Y_hat = @(b,x) b(1).*exp(-b(2)./x);
slope_hat = @(b,x) b(1) .*b(2) .*exp(-b(2)./x) ./ (x.^2);

y_inherent = Y_hat(b_inherent,x_data);

slope_inherent = slope_hat(b_inherent,x_data);

s_inherent = 0.001 * mean(abs(y_inherent));

valid = 0;
valid_slope = 0;

for a= 1:500
    
    b_guess = b_inherent;
    
    Y = normrnd(y_inherent,s_inherent);
        
    SSE_hat = @(b) sum((Y_hat(b,x_data) - Y).^2);
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    Y_fit = Y_hat(b_fit,x_data);
    
    SSE_fit = SSE_hat(b_fit);
    
    sigma2 = SSE_fit / (n - p);
    
    RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));
    
    % This is good for 100 Temperature values between 50 and 500 with 2%
    % error on the mean value
    b0_range = linspace(0.9975,1.0025,100)*b_fit(1);
    b1_range = linspace(0.9975,1.0025,100)*b_fit(2);

    k = 1;

    for i = 1:length(b0_range)
    
        for j=1:length(b1_range)
        
            b_range = [b0_range(i); b1_range(j)];
                    
            SSE_range = SSE_hat(b_range);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_hat(b_range,x_data);
                    
                    slope_ext(k,:) = slope_hat(b_range,x_data);
                
                    k = k+1;
                    
                    if i == 1 || i == length(b0_range) || j == 1 || j == length(b1_range)
                       
                        pause(0)
                        
                    end
                                             
                end
            
        end
    
    end
    
    Y_hi = max(Y_ext);
    Y_lo = min(Y_ext);
    slope_hi = max(slope_ext);
    slope_lo = min(slope_ext);
        
    valid_all = 1;
    valid_slope_all = 1;
        
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < y_inherent(r) || Y_lo(r) > y_inherent(r)
          
           valid_all = 0;
           
       end
       
       if slope_hi(r) < slope_inherent(r) || slope_lo(r) > slope_inherent(r)
           
           valid_slope_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
    end
    
    if valid_slope_all == 1
        
        valid_slope = valid_slope + 1;
        
    end
    
%     figure
%     scatter(b_ext(:,1),b_ext(:,2))

    clear Y_range b_ext Y_ext
    
end

valid/a
valid_slope/a

figure
hold
plot(x_data,y_inherent,'g')
scatter(x_data,Y)
plot(x_data,Y_fit)
plot(x_data,Y_hi,'--')
plot(x_data,Y_lo,'--')
hold


