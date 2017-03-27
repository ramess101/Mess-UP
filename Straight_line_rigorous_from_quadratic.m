clear

N = 1;
x_quad = linspace(-10,10,20)';
x_line = linspace(0.9,1.5,20)';
n_quad = length(x_quad);
n_line = length(x_line);
p_quad = 3;
p_line = 2;

X_quad = ones(n_quad,1);
X_line  = ones(n_line,1);


X_quad = [X_quad, x_quad, x_quad.^2];
X_line_quad = [X_line, x_line, x_line.^2];
X_line = [X_line, x_line];


b_inherent = [2.5; 3.2; 1];

y_inherent_quad = X_quad*b_inherent;
y_inherent_line_quad = X_line_quad * b_inherent;

s_inherent = 0.1;

valid = 0;
valid1 = 0;
valid2 = 0;
valid3 = 0;
valid4 = 0;

for a = 1:100
    
    Y = normrnd(y_inherent_quad,s_inherent);
    
    b_quad = (X_quad'*X_quad)\(X_quad'*Y);
       
    Y_line_quad = X_line_quad*b_quad;
    Y_line_quad = normrnd(Y_line_quad,s_inherent);
    
    b_line = (X_line'*X_line)\(X_line'*Y_line_quad);
    
    Y_hat = X_line * b_line;
    
    SSE_fit = sum((Y_line_quad - Y_hat).^2);
       
    RHS = SSE_fit * (1 + p_quad/(n_quad-p_quad) * (finv(0.95,p_quad,n_quad-p_quad)));
%     RHS = SSE_fit * (1 + p_line/(n_line-p_line) * (finv(0.95,p_line,n_line-p_line)));
    
    % This is good for 20 n_quad and 4 n_line
    b0_range = linspace(0.6,1.4,200)*b_line(1);
    b1_range = linspace(0.6,1.4,200)*b_line(2);

    k = 1;

    for i = 1:length(b0_range)
    
        for j=1:length(b1_range)
        
            b_range = [b0_range(i); b1_range(j)];
        
            SSE_range = sum((Y_line_quad - (X_line*b_range)).^2);
        
                if SSE_range < RHS
                       
                    b_ext(k,:) = b_range;
                
                    k = k+1;
                    
                    if i == 1 || i == length(b0_range) || j == 1 || j == length(b1_range)
                       
                        pause(0)
                        
                    end
                                             
                end
            
        end
    
    end
    
    Y_range = X_line*b_ext';

    Y_hi = max(Y_range')';
    Y_lo = min(Y_range')';
    
    valid_all = 1;
        
    for r = 1:length(Y_hi)
        
       if Y_hi(r) < y_inherent_line_quad(r) || Y_lo(r) > y_inherent_line_quad(r)
          
           valid_all = 0;
           
       end
        
    end
    
    if valid_all == 1
        
        valid = valid + 1;
        
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
    
    clear Y_range b_ext
    
end

valid/a
valid1/a
valid2/a
valid3/a
valid4/a
