clear

x = [0.1 0.1 0.1 0.1 0.1 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.7 0.7 0.7 0.7 0.7]';
n = length(x);

X = ones(n,1);

X = [X, x];

b_inherent = [2.5; 3.2];

y_inherent = X*b_inherent;

s_inherent = 0.1;

t = 2.10092; % Matlab doesn't have a two tail t-statistic

x_avg = mean(x);

x_x2 = (x-x_avg).^2;

trad = s_inherent * t * (1/n + (x_x2/sum(x_x2)));

valid = 0;

for a = 1:100000
    
    Y = normrnd(y_inherent,s_inherent);
    
    b_fit = (X'*X)\(X'*Y);
    
    Y_hat = X*b_fit;
    
    Y_hi = Y_hat + trad;
    Y_lo = Y_hat - trad;
    
%     if Y_hi(1) > y_inherent(1) && Y_lo(1) < y_inherent(1) && Y_hi(6) > y_inherent(6) && Y_lo(6) < y_inherent(6) && Y_hi(11) > y_inherent(11) && Y_lo(11) < y_inherent(11) && Y_hi(16) > y_inherent(16) && Y_lo(16) < y_inherent(16)
%         
%        valid = valid + 1; 
%         
%     end

    if Y_hi(16) > y_inherent(16) && Y_lo(16) < y_inherent(16) 
        
       valid = valid + 1; 
        
    end
    
    
end

valid/a
