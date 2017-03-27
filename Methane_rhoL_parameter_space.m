clear

% From Siepmann's work with GEMC

TC_star = 1.3128;
rhoC_star = 0.316;
PC_star = 0.1274;

T_CH4 = [90.69	94	100	104	110	114	120	124	130	134	140	144	150	154	160	164	170	174	180	184	188]; 
V_CH4 = [28.142	27.866	27.357	27.01	26.478	26.113	25.55	25.163	24.561	24.144	23.491	23.034	22.309	21.794	20.963	20.36	19.354	18.591	17.218	16.043	14.27]; % mol/L

% Using non-weighted regression
% V_fit = [26.523825 26.150069 25.573144 25.17649 24.561024 24.135329 23.47003 23.005948 22.272966  21.755055 20.923163 20.322503 19.327439 18.576918 17.237862 16.085103 14.386061]; % My fit to the data used between 5:21 
% Using relative weighted regression
% V_fit = [26.547061 26.17068 25.589826 25.190561 24.571188 24.142897 23.473721 23.007066  22.270248 21.749799 20.914142 20.311011 19.312331 18.559498 17.217275 16.062872 14.363538];

% These are the values that are reported by McDonald
sig_CH4 = 3.743; %Angstroms
eps_CH4 = 149.1; %Kelvins

T_inherent = linspace(90,190,30);

V_inherent = rhoL_GEMC(T_inherent,eps_CH4,sig_CH4/10);

s_inherent = 0.01 * mean(abs(V_inherent));

valid = 0;

for a= 1:1
    
    b_guess = [eps_CH4,sig_CH4];
    
    % If testing the parameter space code with the inherent values
    
%     T = T_inherent;
%     
%     Y = normrnd(V_inherent,s_inherent);
    
%     If you want to see what the DIPPR uncertainty looks like
    T = T_CH4(5:21);
    Y = V_CH4(5:21); 
    
    [T, I] = sort(T);
    Y = Y(I);
    
%     Y = V_fit;
    
    Y_hat = @(b) rhoL_GEMC(T,b(1),b(2)/10);
    
%     SSE_hat = @(b) sum((Y_hat(b) - Y).^2);
    SSE_hat = @(b) sum(((Y_hat(b) - Y)./Y).^2); % This approach gives a result closer to what is reported
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    Y_fit = rhoL_GEMC(T,b_fit(1),b_fit(2)/10);
    
    SSE_fit = SSE_hat(b_fit);
    
    n = length(T);
    p = 2;
    
    alpha = 0.95;
    
    % Type A analysis
    sigma2 = SSE_fit / (n - p); 
    nu = n-p;
    
    % Type AB analysis
    U_DIPPR = 0.01;
    sigma2 = (U_DIPPR * mean(abs(Y)) / 1.96)^2; 
    nu = 10^6;   
    
    % Type B analysis (relative error)
    Y_hi = Y*(1+U_DIPPR);
    Y_lo = Y*(1-U_DIPPR);
%     
    % Type B analysis (constant absolute error)
%     Y_hi = Y + mean(abs(Y))*U_DIPPR;
%     Y_lo = Y - mean(abs(Y))*U_DIPPR;
%     
    RHS = sigma2 * (n + p * (finv(alpha,p,nu)-1));
       
%     %     This is good for the DIPPR data (5:18) with weighting by relative error
%     eps_range = linspace(0.9975,1.0025,100)*b_fit(1);
%     sig_range = linspace(0.9992,1.0008,100)*b_fit(2);

% %     %     This is good for the DIPPR data (5:21) with weighting by relative error
%     eps_range = linspace(0.9994,1.0006,100)*b_fit(1);
%     sig_range = linspace(0.9995,1.0005,100)*b_fit(2);
    
    %     %     This is good for the DIPPR data (5:21) but at a higher
    %     confidence level for scanning PDF
    eps_range = linspace(0.999,1.001,100)*b_fit(1);
    sig_range = linspace(0.999,1.001,100)*b_fit(2);

        %     This is good for all the DIPPR data with weighting by relative error
%     eps_range = linspace(0.9992,1.0008,100)*b_fit(1);
%     sig_range = linspace(0.9994,1.0006,100)*b_fit(2);

%     %    A more complete scan to make sure there are no local minima
%     eps_range = linspace(0.995,1.005,200)*b_fit(1);
%     sig_range = linspace(0.995,1.005,200)*b_fit(2);

%     %     This is good for the DIPPR data (1:9) with weighting by relative error
%     eps_range = linspace(0.9975,1.0025,100)*b_fit(1);
%     sig_range = linspace(0.9996,1.0004,100)*b_fit(2);

 %    Type B analysis
 % Relative error
 
    eps_range = linspace(0.9975,1.003,400)*b_fit(1);
    sig_range = linspace(0.997,1.004,400)*b_fit(2);
    
     %    Type B analysis
 % Relative error
 % Fit to correlation
%  
%     eps_range = linspace(0.997,1.0032,200)*b_fit(1);
%     sig_range = linspace(0.996,1.004,200)*b_fit(2);
    
 % Absolute error
 
%     eps_range = linspace(0.997,1.004,100)*b_fit(1);
%     sig_range = linspace(0.997,1.004,100)*b_fit(2);
    
  % Absolute error
  % Fit to correlation
 
%     eps_range = linspace(0.997,1.004,200)*b_fit(1);
%     sig_range = linspace(0.9965,1.0035,200)*b_fit(2);
    
    k = 1;
    
    PDF = zeros(length(eps_range),length(sig_range));
    TC_range = PDF;
    rhoC_range = PDF;
    PC_range = PDF;
    eps_sig3_range = PDF;
    PDF_eps = zeros(length(eps_range));
    PDF_sig = zeros(length(sig_range));
    SSE = PDF;

    for i = 1:length(eps_range)
    
        for j=1:length(sig_range)
        
            b_range = [eps_range(i); sig_range(j)];
            
            Y_range = rhoL_GEMC(T,b_range(1),b_range(2)/10);
        
%             SSE_range = SSE_hat(b_range);
%         
%                 if SSE_range < RHS
%                        
%                     b_ext(k,:) = b_range;
%                     
%                     Y_ext(k,:) = Y_range;
%                 
%                     k = k+1;
%                     
%                     if i == 1 || i == length(eps_range) || j == 1 || j == length(sig_range)
%                        
%                         pause(0)
%                         
%                     end
%                                              
%                 end
%                 
%                 PDF(i,j) = fcdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p); % Really should be CDF
%                 PDF_eps(i) = PDF_eps(i) + fpdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
%                 PDF_sig(j) = PDF_sig(j) + fpdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
%                 
%                 SSE(i,j) = SSE_range;

                TC_range(i,j) = eps_range(i) * TC_star; 
                rhoC_range(i,j) = rhoC_star ./ ((sig_range(j)./10).^3) * 1.6605778811;
                PC_range(i,j) = PC_star * eps_range(i) ./ ((sig_range(j)./10).^3) * 0.013806505;
                eps_sig3_range(i,j) = eps_range(i) / (sig_range(j)^3);

        % For Type B analysis
               
            Dev_hi = Y_hi - Y_range;
            Dev_lo = Y_lo - Y_range;
            
%             Max_dev = max(abs(Y-Y_range)); % This approach doesn't let the optimal fit have a CDF of 0
            Max_dev = max(abs(Y_fit-Y_range)); % This way the best fit has a CDF of 0
            
%             t_test = Max_dev/sqrt(sigma2);
                 
            t_test = max(abs(Y_fit-Y_range)./(U_DIPPR*abs(Y_fit)/1.96)); % Here I am saying that it is the largest t_statistic that I want if uncertainty is relative
     
            DOF = 10^6;
                       
            PDF(i,j) = 1 - 2*tcdf(t_test,DOF,'upper'); %Because only returns one-sides t-statistic I needed to be creative
            PDF_eps(i) = PDF_eps(i) + tpdf(t_test,DOF);
            PDF_sig(j) = PDF_sig(j) + tpdf(t_test,DOF);

            if min(Dev_hi) >= 0 && max(Dev_lo) <= 0
                
%                 if eps_range(i) <= b_fit(1) % This includes constraint that fit in TC will not be worse
                                
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_range;
                
                    k = k+1;
                
%                 end

            end     
            
        end
    
    end
    
    Y_hat_hi = max(Y_ext);
    Y_hat_lo = min(Y_ext);
        
%     valid_all = 1;
%         
%     for r = 1:length(Y__hat_hi)
%         
%        if Y_hat_hi(r) < B_inherent(r) || Y_hat_lo(r) > B_inherent(r)
%           
%            valid_all = 0;
%            
%        end
%         
%     end
%     
%     if valid_all == 1
%         
%         valid = valid + 1;
%         
%     end
    
    figure
    hold
    scatter(b_ext(:,2),b_ext(:,1))
    scatter(sig_CH4,eps_CH4,'g')
    scatter(b_fit(2),b_fit(1),'r')
    hold
    
    max(b_ext(:,1))/b_fit(1)
    min(b_ext(:,1))/b_fit(1)
    max(b_ext(:,2))/b_fit(2)
    min(b_ext(:,2))/b_fit(2)
    
%     clear Y_ext b_ext
    
end

PDF_eps = PDF_eps/sum(PDF_eps);
PDF_sig = PDF_sig/sum(PDF_sig);

% valid/a

% file = fopen('B_95.txt','w');
% fwrite(file,valid/a)
% fclose(file);

figure
hold
plot(T_inherent,V_inherent,'g')
scatter(T,Y)
errorbar(T,Y,Y*U_DIPPR) % If relative error
% errorbar(T,Y,Y./Y*mean(abs(Y))*U_DIPPR) % If absolute error
plot(T,Y_fit)
plot(T,Y_hat_hi,'--')
plot(T,Y_hat_lo,'--')
scatter(T_CH4,V_CH4)
hold

figure
contour(sig_range,eps_range,PDF)

figure
plot(eps_range,PDF_eps)

figure
plot(sig_range,PDF_sig)

% histogram_methane
% 
% [eps_counts, eps_centers] = hist(eps_norm,100);
% 
% [eps_low,eps_high,actual_eps] = integrate_histogram(eps_counts,eps_centers,alpha);
% 
% eps_uncertainty = (eps_high - eps_low)/(eps_high+eps_low) * 100;
% 
% [sig_counts, sig_centers] = hist(sig_norm,100);
% 
% [sig_low,sig_high,actual_sig] = integrate_histogram(sig_counts,sig_centers,alpha);
% 
% sig_uncertainty = (sig_high - sig_low)/(sig_high+sig_low) * 100;
% 
% [TC_counts, TC_centers] = hist(TC_norm,100);
% 
% [TC_low,TC_high,actual_TC] = integrate_histogram(TC_counts,TC_centers,alpha);
% 
% TC_uncertainty = (TC_high - TC_low)/(TC_high+TC_low) * 100;
% 
% [rhoC_counts, rhoC_centers] = hist(rhoC_norm,100);
% 
% [rhoC_low,rhoC_high,actual_rhoC] = integrate_histogram(rhoC_counts,rhoC_centers,alpha);
% 
% rhoC_uncertainty = (rhoC_high - rhoC_low)/(rhoC_high+rhoC_low) * 100;
% 
% [PC_counts, PC_centers] = hist(PC_norm,100);
% 
% [PC_low,PC_high,actual_PC] = integrate_histogram(PC_counts,PC_centers,alpha);
% 
% PC_uncertainty = (PC_high - PC_low)/(PC_high+PC_low) * 100;
% 
% % TC_ext = b_ext(:,1) * TC_star;
% % rhoC_ext = rhoC_star ./ ((b_ext(:,2)./10).^3) * 1.6605778811;
% % PC_ext = PC_star * b_ext(:,1) ./ ((b_ext(:,2)./10).^3) * 0.013806505;
% % 
% TC = b_fit(1) * TC_star;
% rhoC = rhoC_star ./ ((b_fit(2)./10).^3) * 1.6605778811;
% PC = PC_star * b_fit(1) ./ ((b_fit(2)./10).^3) * 0.013806505;
% 
% % TC_uncertainty = (max(TC_ext) - min(TC_ext)) / (max(TC_ext) + min(TC_ext)) * 100;
% % rhoC_uncertainty = (max(rhoC_ext) - min(rhoC_ext)) / (max(rhoC_ext) + min(rhoC_ext)) * 100;
% % PC_uncertainty = (max(PC_ext) - min(PC_ext)) / (max(PC_ext) + min(PC_ext)) * 100;

T_plot = linspace(min(T),max(T),1000);
Y_plot_fit = rhoL_GEMC(T_plot,b_fit(1),b_fit(2)/10);
MW=16.0425;
Y_plot_fit = Y_plot_fit * MW / 1000;
Y_plot_low = Y_plot_fit * (1-U_DIPPR);
Y_plot_high = Y_plot_fit * (1+U_DIPPR);

Y_hat_plot_lo = Y_hat_lo * MW / 1000;
Y_hat_plot_hi = Y_hat_hi * MW / 1000;
