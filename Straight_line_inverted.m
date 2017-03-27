clear

N = 1;
x = [0.1*ones(N,1); 0.3*ones(N,1); 0.5*ones(N,1); 0.7*ones(N,1)];
n = length(x);
p = 2;

X = ones(n,1);

X = [X, x];

b_inherent = [2.5; 3.2];

y_inherent = X*b_inherent;

s_inherent = 0.1;
t = 3.182446; % This is using 5 replicates
z = 1.96; % This is at the 95% confidence level
z = 2.49; % This is at the 95^1/4 confidence level

RHS = s_inherent^2 * (n + p * (finv(0.95^(1/4),p,n-p)-1));

b0_range = linspace(0.99,1.01,100)*b_inherent(1);
b1_range = linspace(0.99,1.01,100)*b_inherent(2);

valid = zeros(length(b0_range),length(b1_range));
valid1 = valid;
valid2 = valid;
valid3 = valid;
valid4 = valid;
valid5 = valid;

for i = 1:length(b0_range)
    
   for j=1:length(b1_range)
       
       b_range = [b0_range(i); b1_range(j)];
       
       y_range = X*b_range;

        for a = 1:10000
    
    Y = normrnd(y_range,s_inherent);
       
    SSE_range(a) = sum((Y - y_inherent).^2);
    
    sigma2(a) = SSE_range(a) / (n - p);
    
    sigma(a) = sqrt(sigma2(a));
      
    Y_avg = [mean(Y(1:N)), mean(Y(N+1:2*N)), mean(Y(2*N+1:3*N)), mean(Y(3*N+1:4*N))];
    
    s = [std(Y(1:N)), std(Y(N+1:2*N)), std(Y(2*N+1:3*N)), std(Y(3*N+1:4*N))];
    
    CI = z * s_inherent / sqrt(N);
%     CI = t * s / sqrt(N);

    Y_hi = Y_avg + CI;
    Y_lo = Y_avg - CI;
    
    if Y_hi(1) > y_inherent(1) && Y_lo(1) < y_inherent(1) && Y_hi(2) > y_inherent(N+1) && Y_lo(2) < y_inherent(N+1) && Y_hi(3) > y_inherent(2*N+1) && Y_lo(3) < y_inherent(2*N+1) && Y_hi(4) > y_inherent(3*N+1) && Y_lo(4) < y_inherent(3*N+1)
        
       valid(i,j) = valid(i,j) + 1; 
        
    end

    if Y_hi(1) > y_inherent(1) && Y_lo(1) < y_inherent(1) 
        
       valid1(i,j) = valid1(i,j) + 1; 
        
    end
    
    if Y_hi(2) > y_inherent(N+1) && Y_lo(2) < y_inherent(N+1) 
        
       valid2(i,j) = valid2(i,j) + 1; 
        
    end
    
    if Y_hi(3) > y_inherent(2*N+1) && Y_lo(3) < y_inherent(2*N+1) 
        
       valid3(i,j) = valid3(i,j) + 1; 
        
    end
    
    if Y_hi(4) > y_inherent(3*N+1) && Y_lo(4) < y_inherent(3*N+1) 
        
       valid4(i,j) = valid4(i,j) + 1; 
        
    end
    
    if SSE_range(a) <= RHS
       
        valid5(i,j) = valid5(i,j) + 1;
        
    end

        end
    
    clear y_range
    
    end
    
end

valid/a
% valid1/a
% valid2/a
% valid3/a
% valid4/a

mean(SSE_range)
mean(sigma)

contour(b0_range,b1_range,valid/a)