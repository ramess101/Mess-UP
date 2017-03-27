% Need to run the MC file first to get all the vap_ave and liq_ave

rhoC_range = linspace(0.23,0.243,40);

PDF_rhoC = zeros(length(sigma),length(rhoC_range));

% rhoc_min = 10;
% rhoc_max = 0;

for m = 1:length(sigma)
    
    T = Temp;
    pv = vap_ave(:,m)';
    pL = liq_ave(:,m)';
       
[TC_fit, rhoc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (rhoc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma2 = (SSE_fit/(n-p))^2;

A_range = linspace(0.93,1.07,30)*A_fit;
b_range = linspace(0.99,1.01,30)*b_fit;
TC_range = linspace(0.995,1.005,30)*TC_fit;

for g=1:length(A_range)
    
    for h=1:length(b_range)
    
        for i=1:length(TC_range)
    
            for j=1:length(rhoC_range)
        
        SE = ((y - (rhoC_range(j) + A_range(g)*(TC_range(i)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(i)-T).^beta))./sz).^2;
        
        SSE = sum(SE);
        
        PDF_rhoC(m,j) = PDF_rhoC(m,j) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);

            end 
            
        end
    
    end
    
end

PDF_rhoC(m,:) = PDF_rhoC(m,:)/sum(PDF_rhoC(m,:));

% This helps to find the max and min range you want for the rhoC_range
% Or you can just look at the histograms from before

%     if rhoc_fit > rhoc_max
%     
%         rhoc_max = rhoc_fit;
%     
%     elseif rhoc_fit < rhoc_min
%     
%         rhoc_min = rhoc_fit;
%     
%     end

end

PDF_rhoC_all = sum(PDF_rhoC);

PDF_rhoC_all = PDF_rhoC_all/sum(PDF_rhoC_all);

figure
hold
plot(rhoC_range,PDF_rhoC_all,'k')
hold