clear

% From Siepmann's work with GEMC

TC_star = 1.3128;
rhoC_star = 0.316;
PC_star = 0.1274;

T_CH4 = [110.83	112.43	114.45	116.79	121.25	128.84	136.75	148.28	162.29	178.41	200	202.49	221.1	243.8	250	273.15	273.17	293.15	300	313.15	333.15	353.15	373.15	423.15	473.15	523.15	573.15	623.15	673.15	723.15	773.15	873.15	973.15	1073.15	1173.15	1273.15]; 
B_CH4 = [-0.3301	-0.3199	-0.3078	-0.2955	-0.2745	-0.2443	-0.2189	-0.1877	-0.1584	-0.1322	-0.10542	-0.1034	-0.0858	-0.0703	-0.0662	-0.05443	-0.0537	-0.0457	-0.04302	-0.03826	-0.03184	-0.02623	-0.0213	-0.01125	-3.56E-03	2.46E-03	7.27E-03	0.01115	0.01433	0.01693	0.01909	0.02237	0.02465	0.02625	0.02739	0.02822]; % L/mol 

% These are the data used by Hirschfelder, Curtiss and Bird

T_HCB = [273.15 298.15 323.15 348.15 373.15 398.15 423.15];
B_HCB = [-0.0541162849 -0.043550285 -0.0347894534 -0.0278616881 -0.0216866466 -0.0164366905 -0.0116736764]; 

% These are fit to the values I accepted for my Paper 3 figures

B_fit = [-0.085874 -0.070017 -0.066291 -0.054175 -0.054175 -0.054165 -0.045536 -0.043591 -0.042891 -0.038193 -0.03492 -0.031875 -0.027684 -0.026378 -0.021552 -0.021552 -0.016285 -0.01171];

% These are the values that are reported by McDonald
sig_CH4 = 3.817;
eps_CH4 = 148.2;

T_inherent = linspace(100,1300,1000);

B_inherent = second_virial(T_inherent,eps_CH4,sig_CH4);
% 
% s_inherent = 0.02 * mean(abs(B_inherent));
% 
% valid = 0;

for a= 1:1
    
    b_guess = [eps_CH4,sig_CH4];
    
    % If testing the parameter space code with the inherent values
    
%     T = T_inherent;
%     
%     Y = normrnd(B_inherent,s_inherent);
    
%     If you want to see what the DIPPR uncertainty looks like
    T = T_CH4(13:23);
    Y = B_CH4(13:23); 
    
%     T = [T, T_HCB(7)];
%     Y = [Y, B_HCB(7)];

    T = [T, T_HCB];
    Y = [Y, B_HCB];
    % If you want to use the HCB data
    
%     T = T_HCB;
%     Y = B_HCB;

    [T, I] = sort(T);
    Y = Y(I);
    
    Y = B_fit;
    
    Y_hat = @(b) second_virial(T,b(1),b(2));
    
%     SSE_hat = @(b) sum((Y_hat(b) - Y).^2);
    SSE_hat = @(b) sum(((Y_hat(b) - Y)./Y).^2); % This approach gives a result closer to what is reported
    
    b_fit = fminsearch(SSE_hat,b_guess);
    
    Y_fit = second_virial(T,b_fit(1),b_fit(2));
    
    SSE_fit = SSE_hat(b_fit);
    
    n = length(T);
    p = 2;
    
    alpha = 0.95;
    
    % Type A analysis
    sigma2 = SSE_fit / (n - p); 
    nu = n-p;
    
    % Type AB analysis
    U_DIPPR = 0.05;
    sigma2 = (U_DIPPR * mean(abs(Y)) / 1.96)^2; 
    nu = 10^6;   
    
    % Type B analysis (relative error)
    Y_hi = Y*(1-U_DIPPR); % In this case because Y is negative must define like this
    Y_lo = Y*(1+U_DIPPR);
    
    % Type B analysis (constant absolute error)
%     Y_hi = Y + mean(abs(Y))*U_DIPPR;
%     Y_lo = Y - mean(abs(Y))*U_DIPPR;
%     
    RHS = sigma2 * (n + p * (finv(alpha,p,nu)-1));
       
    %     This is good for the DIPPR data between 13:23 with weighting by relative error
%     eps_range = linspace(0.9935,1.0065,100)*b_fit(1);
%     sig_range = linspace(0.9955,1.0045,100)*b_fit(2);
    
        %     This is good for the DIPPR data between 13:23 with weighting by relative error
%     eps_range = linspace(0.9935,1.0065,100)*b_fit(1);
%     sig_range = linspace(0.994,1.006,100)*b_fit(2);
%     
        %     This is good for the DIPPR data between 13:23 with weighting by relative error
%     eps_range = linspace(0.993,1.007,200)*b_fit(1);
%     sig_range = linspace(0.9935,1.0065,200)*b_fit(2);
%     
            %     This is good for the HCB data with weighting by relative error
%     eps_range = linspace(0.99,1.01,100)*b_fit(1);
%     sig_range = linspace(0.99,1.01,100)*b_fit(2);

        %     This is good for the DIPPR data between 13:23 and the HCB data with weighting by relative error
%     eps_range = linspace(0.996,1.004,100)*b_fit(1);
%     sig_range = linspace(0.996,1.004,100)*b_fit(2);
    
            %  PDF scan, 99.9%   This is good for the DIPPR data between 13:23 and the HCB data with weighting by relative error
%     eps_range = linspace(0.9935,1.0065,200)*b_fit(1);
%     sig_range = linspace(0.9935,1.0065,200)*b_fit(2);

     %     This is good for all the DIPPR data with weighting by relative error
%     eps_range = linspace(0.996,1.004,100)*b_fit(1);
%     sig_range = linspace(0.993,1.007,100)*b_fit(2);

 %     This is good for the DIPPR data between 13:23 and the HCB data with weighting by relative error
 % Type B uncertainty analysis 
 % Constant relative error
%     eps_range = linspace(0.975,1.025,200)*b_fit(1);
%     sig_range = linspace(0.965,1.035,200)*b_fit(2);
%     
     %     This is good for the DIPPR data between 13:23 and the HCB data with weighting by relative error
 % Type B uncertainty analysis 
 % Constant relative error
 % Fit to correlation
    eps_range = linspace(0.975,1.027,100)*b_fit(1);
    sig_range = linspace(0.965,1.035,100)*b_fit(2);
% %     
    %     This is good for the DIPPR data between 13:23 and the HCB data with weighting by relative error
 % Type B uncertainty analysis  
 % Constanta absolute error
%     eps_range = linspace(0.94,1.06,200)*b_fit(1);
%     sig_range = linspace(0.94,1.06,200)*b_fit(2);

    k = 1;
    
    PDF = zeros(length(eps_range),length(sig_range));
    TC_range = PDF;
    rhoC_range = PDF;
    PC_range = PDF;
    PDF_eps = zeros(length(eps_range));
    PDF_sig = zeros(length(sig_range));
    SSE = PDF;

    for i = 1:length(eps_range)
    
        for j=1:length(sig_range)
        
            b_range = [eps_range(i); sig_range(j)];
            
            Y_range = second_virial(T,b_range(1),b_range(2));
            
            TC_range(i,j) = eps_range(i) * TC_star; 
            rhoC_range(i,j) = rhoC_star ./ ((sig_range(j)./10).^3) * 1.6605778811;
            PC_range(i,j) = PC_star * eps_range(i) ./ ((sig_range(j)./10).^3) * 0.013806505;
        
            % For Type A or Type AB analysis
            
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
%                 PDF(i,j) = fcdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
%                 PDF_eps(i) = PDF_eps(i) + fpdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
%                 PDF_sig(j) = PDF_sig(j) + fpdf((SSE_range-SSE_fit)/(sigma2*p),p,n-p);
%                 
%                 SSE(i,j) = SSE_range;

            % For Type B analysis
            
            Dev_hi = Y_hi - Y_range;
            Dev_lo = Y_lo - Y_range;
            
%             Max_dev = max(abs(Y-Y_range)); % This approach doesn't let the optimal fit have a CDF of 0
            Max_dev = max(abs(Y_fit-Y_range)); % This way the best fit has a CDF of 0
            
            t_test = Max_dev/sqrt(sigma2);

%             t_test = max(abs(Y_fit-Y_range)./(U_DIPPR*abs(Y_fit)/1.96)); % Here I am saying that it is the largest t_statistic that I want if uncertainty is relative
                       
            DOF = 10^6;
                       
            PDF(i,j) = 1 - 2*tcdf(t_test,DOF,'upper'); %Because only returns one-sides t-statistic I needed to be creative
            PDF_eps(i) = PDF_eps(i) + tpdf(t_test,DOF);
            PDF_sig(j) = PDF_sig(j) + tpdf(t_test,DOF);
            
            if min(Y_hi-Y_range) >= 0 && max(Y_lo-Y_range) <= 0
                
%                 if eps_range(i) <= b_fit(1) % This includes constraint that fit in TC will not be worse
                
                    b_ext(k,:) = b_range;
                    
                    Y_ext(k,:) = Y_range;
                    
                    Y_ext_plot(k,:) = second_virial(T_inherent,b_range(1),b_range(2));
                
                    k = k+1;
                
%                 end

            end           
                
                
        end
    
    end
    
    Y_hat_hi = max(Y_ext);
    Y_hat_lo = min(Y_ext);
    
    Y_plot_hi = max(Y_ext_plot);
    Y_plot_lo = min(Y_ext_plot);
        
%     valid_all = 1;
%         
%     for r = 1:length(Y_hat_hi)
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
plot(T_inherent,B_inherent,'g')
scatter(T,Y)
% errorbar(T,Y,Y*U_DIPPR) % If relative error
errorbar(T,Y,Y./Y*mean(abs(Y))*U_DIPPR) % If absolute error
plot(T,Y_fit,'r')
plot(T,Y_hat_hi,'--')
plot(T,Y_hat_lo,'--')
plot(T_inherent,Y_plot_hi,'g--')
plot(T_inherent,Y_plot_lo,'g--')
scatter(T_CH4,B_CH4)
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
% % 
% % TC_uncertainty = (max(TC_ext) - min(TC_ext)) / (max(TC_ext) + min(TC_ext)) * 100;
% % rhoC_uncertainty = (max(rhoC_ext) - min(rhoC_ext)) / (max(rhoC_ext) + min(rhoC_ext)) * 100;
% % PC_uncertainty = (max(PC_ext) - min(PC_ext)) / (max(PC_ext) + min(PC_ext)) * 100;

T_plot = linspace(min(T),max(T),1000);
Y_plot_fit = second_virial(T_plot,b_fit(1),b_fit(2));
Y_plot_low = Y_plot_fit * (1+U_DIPPR);
Y_plot_high = Y_plot_fit * (1-U_DIPPR);

% Alternative way for Type A and Type B
%  clear Y_ext, Y_min, Y_max
% for j = 1:length(T)
% for i = 1:length(b_ext)
% Y_ext(i)= second_virial(T(j),b_ext(i,1),b_ext(i,2));
% end
% Y_min(j) = min(Y_ext);
% Y_max(j) = max(Y_ext);
% end
