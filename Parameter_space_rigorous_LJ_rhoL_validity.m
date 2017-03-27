clear

T_data = linspace(90,130,30);

T_LJ = [0.75,0.85,1,1.05,1.1,1.15,1.2];
rhoL_LJ = [0.8195 0.7779 0.701 0.6702, 0.6407, 0.6061, 0.5659];

b_LJ_fit = [0.31443712 0.17405154 1.2996166 0.50021066 0.33325329];

rhoL_LJ_hat = @(T_star) b_LJ_fit(1) + b_LJ_fit(2).*(b_LJ_fit(3) - T_star) + b_LJ_fit(4).*(b_LJ_fit(3) - T_star).^b_LJ_fit(5);

sig_Ar = 0.33873047; % nm
eps_Ar = 115.80313138; % K

T_exp = [87.29 90 90 94.43 99.2 100 100 100 105 105.97 110 110 110  110.98 115 116.81 120 120 120 124.22 130 130 130.03 134.86];
Y_exp = [34.843 34.474 34.483 33.786 33.116 32.834 32.852 32.902 32.031 31.949 31.037 31.063 31.122 30.96 30.162 29.941 29.028 29.048 29.131 28.329 26.669 26.679 26.88 25.253];

Y_hat = @(T, b) rhoL_LJ_hat(T/b(2)) ./ (b(1)^3) * 1.6605778811; % The last value is to convert into mol/L

n = length(T_data); % Make sure this referring to the data that you actually are using, I forgot to change this to T_data from T_exp in one case
p = 2;

b_g = [sig_Ar, eps_Ar];

rhoL_inherent = Y_hat(T_data,b_g);

s_inherent = 0.002 * mean(rhoL_inherent);

valid = 0;

for a = 1:1

rho_data = normrnd(rhoL_inherent, s_inherent);

SSE_hat = @(b) sum((rho_data - Y_hat(T_data,b)).^2);
    
b_fit = fminsearch(SSE_hat,b_g);

Y_fit = Y_hat(T_data,b_fit);
    
SSE_fit = SSE_hat(b_fit);
    
sigma2 = SSE_fit / (n - p);
    
RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));
       
sig_range = linspace(0.998,1.002,80)*b_fit(1);
eps_range = linspace(0.994,1.006,80)*b_fit(2);

    k = 1;
            
for i = 1:length(sig_range)
                
    for j = 1:length(eps_range)
        
            b_range = [sig_range(i); eps_range(j)];
                    
            SSE_range = SSE_hat(b_range);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_hat(T_data,b_range);
                
                    k = k+1;
                    
%                     if SSE_range < SSE_fit
%                         
%                         b_new_best_fit = b_range;
%                         
%                     end
                                                                
                end
                            
    end
    
end
    
    Y_hi = max(Y_ext);
    Y_lo = min(Y_ext);
    
    valid_all = 1;
    
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < rhoL_inherent(r) || Y_lo(r) > rhoL_inherent(r)
          
           valid_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
    end
        
    figure
    plot(b_ext(:,2),b_ext(:,1))
    
%     figure
%     hold
%     scatter(T_data,rhoL_inherent,'g')
% %     scatter(T_data,rho_data)
%     plot(T_data,Y_fit)
%     plot(T_data,Y_hi,'--')
%     plot(T_data,Y_lo,'--')
%     hold

    max(b_ext(:,1))/b_fit(1)
    min(b_ext(:,1))/b_fit(1)
    max(b_ext(:,2))/b_fit(2)
    min(b_ext(:,2))/b_fit(2)
%     
    clear Y_ext b_ext
    
end

valid/a