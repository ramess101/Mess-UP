function [ TC_min,TC_max,pc_min,pc_max ] = rigorous_statistics( T,pv,pL,sv,sL,beta )
%This function requires the best fit parameters as input and it performs a
%rigorous determination of the confidence region.  

% [TC_fit, pc_fit, A_fit, b_fit] = towhee_rectilinear(T,pv,pL,sv,sL,beta);

% sz = ((sv.^2 + sL.^2).^(1/2))/2;
% sy = sz;

[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

% A_fit = 3.8858442*10^-4;
% b_fit = 0.0584264;
% pc_fit = 0.2424774;
% TC_fit = 617.2207397;

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma = SSE_fit/(n-p);

RHS = sigma * (n + p * (finv(0.95,p,n-p)-1)); % No 0.95^p because I want the true 95%

% For the growth part we start with a pretty coarse grid.  Once we have
% grown as much as we want, then we go to a finer grid.

Ab_spacing = 10; % For speed, only 10 values will be used for A and b
pcTC_spacing = 20; % pc and TC are more important
iteration = 0;

% Converged is for the finer grid portion, growing is for the coarse grid
converged = false;
convergence = [false, false, false, false, false, false, false, false];
growing = true;
growth = [true,true,true,true,true,true,true,true];

% Since we are going to grow our confidence region we start with a very
% small region.  These will grow with the coarse grid and then contract
% back with the finer grid.

A_low = 0.99*A_fit;
A_high = 1.01*A_fit;
b_low = 0.995*b_fit;
b_high = 1.005*b_fit;
pc_low = 0.995*pc_fit;
pc_high = 1.005*pc_fit;
% pc_low = 0.2405;
% pc_high = 0.24445;
TC_low = 0.995*TC_fit;
TC_high = 1.005*TC_fit;
% TC_low = 616.75;
% TC_high = 617.735;


while converged ==false
    
A_spacing = (A_high-A_low)/(Ab_spacing-1);
b_spacing = (b_high-b_low)/(Ab_spacing-1);
pc_spacing = (pc_high-pc_low)/(pcTC_spacing-1);
TC_spacing = (TC_high-TC_low)/(pcTC_spacing-1);

A_range = A_low:A_spacing:A_high;
b_range = b_low:b_spacing:b_high;
pc_range = pc_low:pc_spacing:pc_high;
TC_range = TC_low:TC_spacing:TC_high;

Ext = zeros(1,(length(A_range)*length(b_range)*length(pc_range)*length(TC_range)));
A_ext = Ext;
b_ext = Ext;
pc_ext = Ext;
TC_ext = Ext;
clear Ext

s=1;

for g=1:length(A_range)
    
    for h=1:length(b_range)

        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
      
       SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
 
       SSE = sum(SE);
       
       if SSE < RHS 
                 
           A_ext(s) = A_range(g);
           b_ext(s) = b_range(h);
           pc_ext(s) = pc_range(i);
           TC_ext(s) = TC_range(j);
           
           s=s+1;
           
           if SSE < SSE_fit % This is to make sure that we actually have the global minimum. If not, rerun the analysis with the new optimum (after plugging in as a new Mathcad guess)
               
              new_best_fit = [A_range(g) b_range(h) pc_range(i) TC_range(j)]; 
              SSE_fit = SSE;
              
           end
           
       end
       
            end
            
       end
       
   end
    
end

% This eliminates the superfluous zero elements
A_ext = A_ext(:,1:(s-1));
b_ext = b_ext(:,1:(s-1));
pc_ext = pc_ext(:,1:(s-1));
TC_ext = TC_ext(:,1:(s-1));

A_low_temp = min(A_ext);
A_high_temp = max(A_ext);
b_low_temp = min(b_ext);
b_high_temp = max(b_ext);
pc_low_temp = min(pc_ext);
pc_high_temp = max(pc_ext);
TC_low_temp = min(TC_ext);
TC_high_temp = max(TC_ext);

if growing % The first part is to grow the regions
    
    if isempty(A_low_temp)
        
        factor = 10^(-(6+iteration));
        low_factor = 1-factor;
        high_factor = 1+factor;
        
        A_low = low_factor*A_fit;
        A_high = high_factor*A_fit;
        b_low = low_factor*b_fit;
        b_high = high_factor*b_fit;
        pc_low = low_factor*pc_fit;
        pc_high = high_factor*pc_fit;
        TC_low = low_factor*TC_fit;
        TC_high = high_factor*TC_fit;
        
        
    else
    
    if A_low_temp == A_low  
    
    A_low = A_low - 2*A_spacing;
    growth(1) = true;
    
    else
    
    growth(1) = false;
    
    end

    if A_high_temp == A_high  
    
    A_high = A_high + 2*A_spacing;
    growth(2) = true;
      
    else
    
    growth(2) = false;
    
    end
    
    if b_low_temp == b_low  
    
    b_low = b_low - 2*b_spacing;
    growth(3) = true;
    
    else
    
    growth(3) = false;
    
    end

    if b_high_temp == b_high  
    
    b_high = b_high + 2*b_spacing;
    growth(4) = true;
    
    else
    
    growth(4) = false;
    
    end

    if pc_low_temp == pc_low  
    
    pc_low = pc_low - 2*pc_spacing;
    growth(5) = true;
    
    else
    
    growth(5) = false;
    
    end

    if pc_high_temp == pc_high  
    
    pc_high = pc_high + 2*pc_spacing;
    growth(6) = true;
    
    else
    
    growth(6) = false;
    
    end

    if TC_low_temp == TC_low  
    
    TC_low = TC_low - 2*TC_spacing;
    growth(7) = true;
    
    else
    
    growth(7) = false;
    
    end

    if TC_high_temp == TC_high  
    
    TC_high = TC_high + 2*TC_spacing;
    growth(8) = true;
    
    else
    
    growth(8) = false;
    
    end
    
    end
    
    if growth(:)==false % Once all the regions incapsiluate the extrema we contract back down with a more refined grid
        
    growing=false;
    Ab_spacing = 20; % For speed, only 20 values will be used for A and b
    pcTC_spacing = 40; % pc and TC are more important
    
    end
    
else

if isempty(A_low_temp) % If one is empty they are all empty (could also use if s==1)
    
    A_low = A_low + A_spacing/2;
    A_high = A_high - A_spacing/2;
    b_low = b_low + b_spacing/2;
    b_high = b_high - b_spacing/2;
    pc_low = pc_low + pc_spacing/2;
    pc_high = pc_high - pc_spacing/2;
    TC_low = TC_low + TC_spacing/2;
    TC_high = TC_high - TC_spacing/2;
    
else

    if A_low_temp == A_low  
    
   convergence(1) = true;  
    
    else
    
   A_low = A_low_temp - A_spacing/2;
    
    end

    if A_high_temp == A_high  
    
   convergence(2) = true;
      
    else
    
   A_high = A_high_temp + A_spacing/2;
    
    end
    
    if b_low_temp == b_low  
    
   convergence(3) = true;
    
    else
    
    b_low = b_low_temp - A_spacing/2;
    
    end

    if b_high_temp == b_high  
    
    convergence(4) = true;
    
    else
    
    b_high = b_high_temp + A_spacing/2;
    
    end

    if pc_low_temp == pc_low  
    
    convergence(5) = true;
    
    else
    
    pc_low = pc_low_temp - pc_spacing/2;
    
    end

    if pc_high_temp == pc_high  
    
    convergence(6) = true;
    
    else
    
    pc_high = pc_high_temp + pc_spacing/2;
    
    end

    if TC_low_temp == TC_low  
    
    convergence(7) = true;
    
    else
    
    TC_low = TC_low_temp - TC_spacing/2;
    
    end

    if TC_high_temp == TC_high  
    
    convergence(8) = true;
    
    else
    
    TC_high = TC_high_temp + TC_spacing/2;
    
    end

end

end

if convergence(:)==true % If all have converged we are done
        
    converged=true;
    
end
    
iteration = iteration + 1;

if iteration == 30 % To avoid infinite loops
    
    converged = true;
    
end

end

subplot(3,2,1)
plot(A_ext,b_ext)
subplot(3,2,2)
plot(pc_ext,b_ext)
subplot(3,2,3)
plot(A_ext,TC_ext)
subplot(3,2,4)
plot(pc_ext,TC_ext)
subplot(3,2,5)
plot(A_ext,pc_ext)
subplot(3,2,6)
plot(b_ext,TC_ext)

% % % To plot just the pc and TC region
% 
% pc_scan=min(pc_range):pc_spacing:max(pc_range);
% 
% k=1;
% 
% for h=1:length(pc_scan)
%     
% j=1;
% 
%     for i=1:length(pc_ext)
%     
% 	if pc_ext(i) == pc_scan(h)
% 
% %         if (pc_ext(i) - pc_scan(h)) < 0.00000001
%     
%     TC_scan(j) = TC_ext(i);
% 
%     j=j+1;
% 
%         end
% 
%     end
% 
% % I did it this way because, for some reason, not every pc in the range of 
% % pc_scan actually has an accepted point.  Somehow I missed the A,B and TC 
% % combinations that were required.
% 
%     if j>1
%        
%     TC_upper(k) = max(TC_scan);
%     TC_lower(k) = min(TC_scan);
%     pc_plot(k) = pc_scan(h);
%     
%     k=k+1;
%     
%     end
%     
% TC_scan=TC_fit;
% 
% end
% 
% % hold
% % 
% % plot(pc_ext,TC_ext)
% % plot(pc_plot,TC_upper);
% % plot(pc_plot,TC_lower);
% % 
% % hold

% TC_lower=TC_lower';
% TC_upper=TC_upper';
% pc_plot=pc_plot';

if iteration < 30
    
    TC_min = min(TC_ext);
    TC_max = max(TC_ext);
    pc_min = min(pc_ext);
    pc_max = max(pc_ext);
else
    TC_min = 0;
    TC_max = 0;
    pc_min = 0;
    pc_max = 0;

end

A_low/A_fit
b_low/b_fit
pc_low/pc_fit
TC_low/TC_fit




