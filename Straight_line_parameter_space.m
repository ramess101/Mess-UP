clear

x = [0.1 0.1 0.1 0.1 0.1 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.7 0.7 0.7 0.7 0.7];
% First set of data were from a really good fit, the second set are from a
% poor fit (i.e. large uncertainty)
% y = [2.845699292 2.862378705 2.849124137 2.720128448 2.870524619 3.650652653 3.315452657 3.47320683 3.396843626 3.371151831 4.189169578 4.250385467 4.162804061 4.047507068 4.062835584 4.834365795 4.690618039 4.92142106 4.519304773 4.771340383];
y = [2.785089681 2.865972683 2.988803154 2.795095154 2.841715542 3.494601775 3.409683474 3.459671508 3.334718858 3.611284788 4.126666899  4.10783323 4.204684124 4.220587919 4.031849002 4.868704945 4.919069646 4.791730986 4.616357024 4.917052726]; 

n = length(x);

X = ones(n,1);

X = [X, x'];

Y = y';

b_fit = (X'*X)\(X'*Y);

SSE_fit = sum((Y - (X*b_fit)).^2);

p = 2;

sigma2 = SSE_fit / (n - p);

% sigma2 = 0.1^2;

RHS = sigma2 * (n + p * (finv(0.95,p,n-p)-1));

% RHS = SSE_fit*(1+(p/(n-p))*finv(0.95,p,n-p));

b0_range = linspace(0.94,1.06,100)*b_fit(1);
b1_range = linspace(0.905,1.095,100)*b_fit(2);

k = 1;

for i = 1:length(b0_range)
    
    for j=1:length(b1_range)
        
        b_range = [b0_range(i); b1_range(j)];
        
        SSE_range = sum((Y - (X*b_range)).^2);
        
            if SSE_range < RHS
                       
                b_ext(k,:) = b_range;
                
                k = k+1;
                                             
            end
            
    end
    
end

b0_max = max(b_ext(:,1))/b_fit(1)
b1_max = max(b_ext(:,2))/b_fit(2)
b0_min = min(b_ext(:,1))/b_fit(1)
b1_min = min(b_ext(:,2))/b_fit(2)

Y_range = X*b_ext';

Y_max = max(Y_range');
Y_min = min(Y_range');

hold
scatter(2.5,3.2,'r')
scatter(b_ext(:,1),b_ext(:,2),'b')
hold

% hold
% scatter(x,y)
% plot(x,Y_max)
% plot(x,Y_min)
% hold