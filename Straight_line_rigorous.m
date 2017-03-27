clear

N = 5;
x = [0.1*ones(N,1); 0.3*ones(N,1); 0.5*ones(N,1); 0.7*ones(N,1)];
n = length(x);
p = 2;

X = ones(n,1);

X = [X, x];
X_quad = [X, x.^2];

b_line = [2.5; 3.2]; % Truely linear
b_quad = [2.5; 3.2; 0.01]; % Not actually linear

y_line = X*b_line;
y_quad = X_quad*b_quad;  % For a not truly linear set of data

y_inherent = y_line;
slope_inherent = b_line(2);

s_inherent = 0.1;
alpha = 0.95;

valid = 0;
valid1 = 0;
valid2 = 0;
valid3 = 0;
valid4 = 0;
valid_slope = 0;

for a = 1:1000
    
    Y = normrnd(y_inherent,s_inherent);
    
    b_fit = (X'*X)\(X'*Y);
    
    Y_hat = X*b_fit;
    
    SSE_fit = sum((Y - Y_hat).^2);
    
    sigma2 = SSE_fit / (n - p);
    
    RHS = sigma2 * (n + p * (finv(alpha,p,n-p)-1));
    
    % This is good for 5 replicates
    b0_range = linspace(0.94,1.06,60)*b_fit(1);
    b1_range = linspace(0.905,1.095,60)*b_fit(2);

    % This is good for 4 replicates
%     b0_range = linspace(0.9,1.1,100)*b_fit(1);
%     b1_range = linspace(0.85,1.15,100)*b_fit(2);
    
    % This is good for 1 replicate
%     b0_range = linspace(-0.1,2,200)*b_fit(1);
%     b1_range = linspace(-0.1,2,200)*b_fit(2);

    % This is good for 1000 replicates
%     b0_range = linspace(0.98,1.02,100)*b_fit(1);
%     b1_range = linspace(0.98,1.02,100)*b_fit(2);
%     
    % This is good for 1000 replicates but with inherently quadratic data
%     b0_range = linspace(0.99,1.01,100)*b_fit(1);
%     b1_range = linspace(0.99,1.01,100)*b_fit(2);

    k = 1;

    for i = 1:length(b0_range)
    
        for j=1:length(b1_range)
        
            b_range = [b0_range(i); b1_range(j)];
        
            SSE_range = sum((Y - (X*b_range)).^2);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                
                    k = k+1;
                    
                    if i == 1 || i == length(b0_range) || j == 1 || j == length(b1_range)
                       
                        pause(0)
                        
                    end
                                             
                end
                
                PDF(i,j) = fcdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
            
        end
    
    end
    
    Y_range = X*b_ext';

    Y_hi = max(Y_range');
    Y_lo = min(Y_range');
    
    valid_all = 1;
        
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < y_inherent(r) || Y_lo(r) > y_inherent(r)
          
           valid_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
    end
    
    % The traditional approach
%     slope_hi = max(b_ext(:,2));
%     slope_lo = min(b_ext(:,2));

%     % The rigorous approach
    for i = 1:1000
        
%         r=normcdf(normrnd(0,1),0,1)*ones(length(b0_range),length(b1_range));
    
        r=fcdf(frnd(p,n-p),p,n-p)*ones(length(b0_range),length(b1_range));
    
        dif = (r-PDF).^2;
    
        [dif_min, I] = min(dif);
    
        [dif_min, J] = min(dif_min);
    
        b0_norm(i) = b0_range(I(J));
        
        b1_norm(i) = b1_range(J);
    
    end
    
    [b1_counts, b1_centers] = hist(b1_norm,100);
    
    [slope_lo, slope_hi] = integrate_histogram(b1_counts,b1_centers,alpha);
    
    if slope_hi >= slope_inherent && slope_lo <= slope_inherent
        
        valid_slope = valid_slope + 1;
        
    end
    
%     if Y_hi(1) > y_inherent(1) && Y_lo(1) < y_inherent(1) && Y_hi(N+1) > y_inherent(N+1) && Y_lo(N+1) < y_inherent(N+1) && Y_hi(2*N+1) > y_inherent(2*N+1) && Y_lo(2*N+1) < y_inherent(2*N+1) && Y_hi(3*N+1) > y_inherent(3*N+1) && Y_lo(3*N+1) < y_inherent(3*N+1)
%         
%        valid = valid + 1; 
%         
%     end
% 
%     if Y_hi(1) > y_inherent(1) && Y_lo(1) < y_inherent(1) 
%         
%        valid1 = valid1 + 1; 
%         
%     end
%     
%     if Y_hi(N+1) > y_inherent(N+1) && Y_lo(N+1) < y_inherent(N+1) 
%         
%        valid2 = valid2 + 1; 
%         
%     end
%     
%     if Y_hi(2*N+1) > y_inherent(2*N+1) && Y_lo(2*N+1) < y_inherent(2*N+1) 
%         
%        valid3 = valid3 + 1; 
%         
%     end
%     
%     if Y_hi(3*N+1) > y_inherent(3*N+1) && Y_lo(3*N+1) < y_inherent(3*N+1) 
%         
%        valid4 = valid4 + 1; 
%         
%     end

%     hold
%     scatter(x,y_quad)
%     plot(x,Y_hi,'--')
%     plot(x,Y_lo,'--')
%     hold
% 
%     figure
%     scatter(b_ext(:,1),b_ext(:,2))
% 
%     max(b_ext(:,1))/b_fit(1)
%     min(b_ext(:,1))/b_fit(1)
%     max(b_ext(:,2))/b_fit(2)
%     min(b_ext(:,2))/b_fit(2)
    
    clear Y_range b_ext
    
end

valid/a
% valid1/a
% valid2/a
% valid3/a
% valid4/a
valid_slope/a
