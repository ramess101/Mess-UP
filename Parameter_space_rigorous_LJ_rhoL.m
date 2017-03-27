clear

T_LJ = [0.75,0.85,1,1.05,1.1,1.15,1.2];
rhoL_LJ = [0.8195 0.7779 0.701 0.6702, 0.6407, 0.6061, 0.5659];

b_LJ_fit = [0.31443712 0.17405154 1.2996166 0.50021066 0.33325329];

rhoL_LJ_hat = @(T_star) b_LJ_fit(1) + b_LJ_fit(2).*(b_LJ_fit(3) - T_star) + b_LJ_fit(4).*(b_LJ_fit(3) - T_star).^b_LJ_fit(5);

sig_Ar = 0.33873047; % nm
eps_Ar = 115.80313138; % K

T_exp = [87.29 90 90 94.43 99.2 100 100 100 105 105.97 110 110 110  110.98 115 116.81 120 120 120 124.22 130 130 130.03 134.86];
Y_exp = [34.843 34.474 34.483 33.786 33.116 32.834 32.852 32.902 32.031 31.949 31.037 31.063 31.122 30.96 30.162 29.941 29.028 29.048 29.131 28.329 26.669 26.679 26.88 25.253];

DIPPR = [34.2901, 32.2206, 28.5972, 27.1579, 25.5144, 23.5257, 20.7578]; 

Y_hat = @(sig,rhoL_LJ) rhoL_LJ ./ (sig^3) * 1.6605778811; % The last value is to convert this into mol/L
X_hat = @(eps,T_LJ) T_LJ .* eps;

rhoL_Ar = Y_hat(sig_Ar,rhoL_LJ);
T_Ar = X_hat(eps_Ar,T_LJ);

T_temp = T_Ar;
rhoL_temp = rhoL_Ar;

Y_hat = @(T, sig, eps) rhoL_LJ_hat(T/eps) ./ (sig^3) * 1.6605778811; % The last value is to convert into mol/L

SSE = @(sig,eps) sum((Y_exp - Y_hat(T_exp,sig,eps)).^2);

N = 5;

for i = 1:N-1
   
    T_temp = [T_temp, T_Ar];
    rhoL_temp = [rhoL_temp, rhoL_Ar];
    
end

T_Ar = T_temp;
rhoL_Ar = rhoL_temp;

n = length(T_Ar);
p = 2;

y_inherent = rhoL_Ar;
x_inherent = T_Ar;

sy_inherent = 0.827; % Constant uncertainty
sx_inherent = 0.1; 

valid = 0;

for a= 1:1
    
    b_guess = [eps_Ar,sig_Ar];
    
    Y = normrnd(y_inherent,sy_inherent);
    X = normrnd(x_inherent,sx_inherent);
        
    SSE_hat = @(b) sum((Y_hat(b(2),Y) - Y).^2 + (X_hat(b(1),X) - X).^2);
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    X_fit = X_hat(b_fit(1),X);
    
    Y_fit = Y_hat(b_fit(2),Y);
    
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
                    
                    Y_ext(k,:) = Y_hat(b_range,x_inherent);
                
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
plot(x_inherent,y_inherent,'g')
scatter(x_inherent,Y)
plot(x_inherent,Y_fit)
plot(x_inherent,Y_hi,'--')
plot(x_inherent,Y_lo,'--')
hold


