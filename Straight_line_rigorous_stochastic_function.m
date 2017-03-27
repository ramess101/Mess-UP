clear

N = 12;
x_exp = linspace(0.5,2,12)';
x_fun = linspace(0.5,2,6)';
x_fun_temp = x_fun;

for i = 1:N-1
   
    x_fun_temp = [x_fun_temp; x_fun];
    
end

x_fun = x_fun_temp;
clear x_fun_temp

n_exp = length(x_exp);
n_fun = length(x_fun);
n_x = n_fun/N;
p = 2;

X_fun  = ones(n_fun,1);
X_exp = ones(n_exp,1);

X_exp = [X_exp, x_exp];
X_fun = [X_fun, x_fun];

b_inherent = [2.5; 3.2];

y_inherent_exp = X_exp*b_inherent;

s_exp = 0.001;
s_fun = 0.01;

per_exp = s_exp / mean(y_inherent_exp) * 100;
per_fun = s_fun / mean(y_inherent_exp) * 100;

valid = 0;
valid_alt = 0;
SSE_fit_exp_ave = 0;
SSE_fit_ave = 0;
SSE_min_ave = 0;
RHS_ave = 0;
RMS_fit_ave = 0;
RMS_fit_exp_ave = 0;

cycles = 1;

Y_exp = normrnd(y_inherent_exp,s_exp); %Not sure if this should be included in the loop

    b_exact_fit = (X_exp'*X_exp)\(X_exp'*Y_exp);
    
    Y_exact_fit = X_exp*b_exact_fit;
    
    SSE_exact_fit = sum((Y_exact_fit - Y_exp).^2);
    
    RHS_exact_fit = SSE_exact_fit * (1+p/(n_exp-p)*(finv(0.95,p,n_exp-p)));

for a = 1:cycles
    
    b0_range = linspace(0.998,1.002,100)*b_inherent(1);
    b1_range = linspace(0.998,1.002,100)*b_inherent(2);
    
    SSE_fit = zeros(length(b0_range),length(b1_range));
    SSE_fit_exp = SSE_fit;
    RMS_fit = SSE_fit;
    RMS_fit_exp = SSE_fit;
    
    m = 1;
    
    for i = 1:length(b0_range)
        
        for j = 1:length(b1_range)
            
        b_range = [b0_range(i); b1_range(j)];
        
        y_inherent_fun = X_fun*b_range;
        
        Y_fun = normrnd(y_inherent_fun,s_fun);

        b_fit = (X_fun'*X_fun)\(X_fun'*Y_fun);
        
        Y_fun_fit = X_fun*b_fit;
        
        Y_fun_exp = X_exp*b_fit;
        
        Y_fun_ave = zeros(n_x,1);
        
        for q = 1:(n_x)
            
            for t = 1:N
               
                Y_fun_ave(q) = Y_fun_ave(q) + Y_fun(q + n_x*(t-1));
                
            end
            
        end
        
        Y_fun_ave = Y_fun_ave / N;
        
        SSE_fun_ave = sum((Y_fun_fit(1:n_x) - Y_fun_ave).^2); % Just relative to the average
        
        SSE_fun = sum((Y_fun_fit - Y_fun).^2); % The difference between the fit to the stochastic function results and the stochastic results themselves
        
        SSE_fun_exp = sum((Y_fun_exp - Y_exp).^2); % The difference between the fit to the stochastic function results and the experimental data
        
        SSE_fit(i,j) = SSE_fun_ave;
        
        RMS_fit(i,j) = sqrt(SSE_fun/n_fun);
        
        SSE_fit_exp(i,j) = SSE_fun_exp;
        
        RMS_fit_exp(i,j) = sqrt(SSE_fun_exp/n_exp);
        
        % Here I am finding what the real 95% confidence region looks like
        
        Y_exact_fun_exp = X_exp*b_range; % The exact function evaluated at the exp
        
        SSE_exact_fun_exp(i,j) = sum((Y_exact_fun_exp - Y_exp).^2); % The difference between the exact function and the experimental
           
            if a == 1 && SSE_exact_fun_exp(i,j) < RHS_exact_fit
                
                b_ext_exact(m,:) = b_range;
                
                m = m+1;
                
            end        
        
        end
            
    end
    
       % This smoothes the results with a quadratic
       
       SSE_smooth = SSE_fit_exp;
       
       % Loop until converges
       
       for h=1:5
       
       SSE_fit_exp = SSE_smooth; 
           
   for j=1:length(b1_range)
       
           x = b0_range';
           X = ones(length(x),1);
           X = [X,x,x.^2];
           Y = SSE_fit_exp(:,j);
           b = (X'*X)\(X'*Y);
           SSE_smooth(:,j) = X*b;
       
   end
   
   for i = 1:length(b0_range)

           x = b1_range';
           X = ones(length(x),1);
           X = [X,x,x.^2];
           Y = SSE_fit_exp(i,:)';
           b = (X'*X)\(X'*Y);
           SSE_smooth(i,:) = X*b;
       
   end
   
       end
   
    RMS_smooth = sqrt(SSE_smooth/n_exp);
    
    SSE_all = SSE_fit + SSE_fit_exp;
               
    [SSE_min, I] = min(SSE_smooth);
    
    [SSE_min, J] = min(SSE_min);
  
    b_opt = [b0_range(I(J)), b1_range(J)];
    
    RHS = SSE_min * (1 + p/(n_exp-p) * (finv(0.95,p,n_exp-p))); % If I add n_exp to n_x I seem to get a good result. Don't know why
    
    sigma2 = SSE_min/(n_exp-p);
    
    sigma2 = 0.1^2 + 0.1^2;
    
    sigma2 = 0.1^2;
    
%     RHS = sigma2 * (n_exp + p*(finv(0.95,p,n_exp-p) - 1));
    
    k = 1;
    
    if SSE_fit_exp(50,50) <= RHS
        
        valid = valid+1;
        
    end
    
    for i = 1:length(b0_range)
        
        for j = 1:length(b1_range)
            
            b_range = [b0_range(i); b1_range(j)];
            
            if SSE_smooth(i,j) <= RHS
               
                b_ext(k,:) = b_range;
                
                k = k+1;
                
                    if i == 1 || i == length(b0_range) || j == 1 || j == length(b1_range)
                       
                        pause(0)
                        
                    end
                    
%                     if i == 50 && j == 50
%                         
%                         valid = valid + 1;
%                         
%                     end
                
            end
            
        end
        
    end
    
    hold
    scatter(b_ext(:,1),b_ext(:,2),'r')
    hold

Y_ext = X_exp * b_ext';

Y_hi = max(Y_ext')';
Y_lo = min(Y_ext')';

valid_all = 1;

    for r = 1:length(Y_hi)
        
       if Y_hi(r) < y_inherent_exp(r) || Y_lo(r) > y_inherent_exp(r)
          
           valid_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid_alt = valid_alt + 1;
        
    end

RMS_fit_ave = RMS_fit_ave + RMS_fit;
RMS_fit_exp_ave = RMS_fit_exp_ave + RMS_fit_exp;
RHS_ave = RHS_ave + RHS;
SSE_min_ave = SSE_min_ave + SSE_min;
SSE_fit_ave = SSE_fit_ave + SSE_fit;
SSE_fit_exp_ave = SSE_fit_exp_ave + SSE_fit_exp;
     
end

   hold
   scatter(b_ext_exact(:,1),b_ext_exact(:,2),'bs')
   hold
      
   m = 1;
     
   for i = 1:length(b0_range)
       
       for j=1:length(b1_range)
           
                if SSE_exact_fun_exp(i,j) < RHS_ave
                    
                b_range = [b0_range(i); b1_range(j)];
                
                b_ext_exact(m,:) = b_range;
                
                m = m+1;
                
                end   
           
       end
       
   end
   
%    hold
%    scatter(b_ext_exact(:,1),b_ext_exact(:,2),'g*')
%    hold

valid/a
valid_alt/a

RMS_fit_ave = RMS_fit_ave / a;
RMS_fit_exp_ave = RMS_fit_exp_ave / a;
RHS_ave = RHS_ave / a;
RMS_RHS_ave = sqrt(RHS_ave/n_exp);
SSE_min_ave = SSE_min_ave /a;
SSE_fit_ave= SSE_fit_ave / a;
SSE_fit_exp_ave = SSE_fit_exp_ave / a;

mean(mean(RMS_fit_ave))
RMS_fit_exp_ave(50,50)

s_exp/sqrt(n_exp)
s_fun/sqrt(N)

(max(b_ext(:,1)-min(b_ext(:,1))))/max(b_ext_exact(:,1)-min(b_ext_exact(:,1)))
(max(b_ext(:,2)-min(b_ext(:,2))))/max(b_ext_exact(:,2)-min(b_ext_exact(:,2)))

% figure
% contour(b0_range,b1_range,SSE_fit_ave)
% 
% figure
% contour(b0_range,b1_range,RMS_fit_exp_ave)