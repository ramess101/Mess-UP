clear

T_LJ = [0.75,0.85,1,1.05,1.1,1.15,1.2];
rhoL_LJ = [0.8195 0.7779 0.701 0.6702, 0.6407, 0.6061, 0.5659];

b_LJ_fit = [0.31443712 0.17405154 1.2996166 0.50021066 0.33325329];

rhoL_LJ_hat = @(T_star) b_LJ_fit(1) + b_LJ_fit(2).*(b_LJ_fit(3) - T_star) + b_LJ_fit(4).*(b_LJ_fit(3) - T_star).^b_LJ_fit(5);

sig_Ar = 0.33873047; % nm
eps_Ar = 115.80313138; % K

T_exp = [87.29 90 90 94.43 99.2 100 100 100 105 105.97 110 110 110  110.98 115 116.81 120 120 120 124.22 130 130 130.03 134.86];
Y_exp = [34.843 34.474 34.483 33.786 33.116 32.834 32.852 32.902 32.031 31.949 31.037 31.063 31.122 30.96 30.162 29.941 29.028 29.048 29.131 28.329 26.669 26.679 26.88 25.253];

Y_hat = @(T, sig, eps) rhoL_LJ_hat(T/eps) ./ (sig^3) * 1.6605778811; % The last value is to convert into mol/L

SSE_hat = @(sig,eps) sum((Y_exp - Y_hat(T_exp,sig,eps)).^2);

n = length(T_exp);
p = 2;
    
SSE_fit = SSE_hat(sig_Ar,eps_Ar);
    
sigma2 = SSE_fit / (n - p);
    
RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));
       
sig_range = linspace(0.9985,1.0015,200)*sig_Ar;
eps_range = linspace(0.9946,1.0056,200)*eps_Ar;

    k = 1;
            
for i = 1:length(sig_range)
                
    for j = 1:length(eps_range)
        
            b_range = [sig_range(i); eps_range(j)];
                    
            SSE_range = SSE_hat(b_range(1),b_range(2));
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_hat(T_exp,b_range(1),b_range(2));
                
                    k = k+1;
                    
                    if SSE_range < SSE_fit
                        
                        b_new_best_fit = b_range;
                        
                    end
                                                                
                end
                            
    end
    
end
    
    Y_hi = max(Y_ext);
    Y_lo = min(Y_ext);
        
    figure
    plot(b_ext(:,2),b_ext(:,1))

    max(b_ext(:,1))/sig_Ar
    min(b_ext(:,1))/sig_Ar
    max(b_ext(:,2))/eps_Ar
    min(b_ext(:,2))/eps_Ar

% figure
% hold
% plot(x_inherent,y_inherent,'g')
% scatter(x_inherent,Y)
% plot(x_inherent,Y_fit)
% plot(x_inherent,Y_hi,'--')
% plot(x_inherent,Y_lo,'--')
% hold


