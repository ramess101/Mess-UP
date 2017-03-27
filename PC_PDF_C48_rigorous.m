% Need to run the MC file first to get all the vap_ave and liq_ave

MW = 675.2878;

% PC_range = linspace(0.3,0.6,100);
% ZC_range = linspace(0.165,0.23,100);

PC_range = linspace(0.4,0.46,100);
ZC_range = linspace(0.183,0.189,100);

PDF_PC = zeros(length(sigma),length(PC_range));
PDF_ZC = zeros(length(sigma),length(ZC_range));

for m = 1:length(sigma)
    
    T = Temp;
    pv = vap_ave(:,m)';
    pL = liq_ave(:,m)';
       
[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

n = 2*length(T);

TR = T(1)/TC_fit;

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma2 = (SSE_fit/(n-p))^2;

% A_range = linspace(0.75,1.25,30)*A_fit;
% b_range = linspace(0.96,1.04,30)*b_fit;
% pc_range = linspace(0.85,1.15,30)*pc_fit;
% TC_range = linspace(0.975,1.025,30)*TC_fit;

A_range = linspace(0.947,1.053,40)*A_fit;
b_range = linspace(0.992,1.008,40)*b_fit;
pc_range = linspace(0.974,1.026,40)*pc_fit;
TC_range = linspace(0.995,1.005,40)*TC_fit;

% s = 1;

% PDF_c = ones(length(A_range)*length(b_range)*length(pc_range)*length(TC_range),1);
% PDF_p = 0 * PDF_c;
% PC = 0 * PDF_p;
% ZC = PC;

for g=1:length(A_range)
    
    for h=1:length(b_range)
    
        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
        
        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
        
        SSE = sum(SE);

        PC = PC_from_Rackett_D(TC_range(j),pc_range(i),MW/pc_range(i),A_range(g),b_range(h),TR,2/7);
       
        ZC = PC * MW / (TC_range(j) * pc_range(i) * 8.314472);
        
        % Determine which PC and ZC bin
        
        dif_PC = abs(PC-PC_range);
        
        [dif_min_PC,k] = min(dif_PC);
        
        dif_ZC = abs(ZC-ZC_range);
        
        [dif_min_ZC,kk] = min(dif_ZC);
               
        PDF_PC(m,k) = PDF_PC(m,k) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
        
        PDF_ZC(m,kk) = PDF_ZC(m,kk) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
        
        % Must make PC and ZC (s) to use this to check
        
%         PDF_p(s) = fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
%         
%         PDF_c(s) = fcdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
%         
%         s = s+1;
                        
            end 
            
        end
    
    end
    
end

PDF_PC(m,:) = PDF_PC(m,:)/sum(PDF_PC(m,:));

PDF_ZC(m,:) = PDF_ZC(m,:)/sum(PDF_ZC(m,:));

end

PDF_PC_all = sum(PDF_PC);

PDF_PC_all = PDF_PC_all/sum(PDF_PC_all);

figure
hold
plot(PC_range,PDF_PC_all,'k')
hold

PDF_ZC_all = sum(PDF_ZC);

PDF_ZC_all = PDF_ZC_all/sum(PDF_ZC_all);

figure
hold
plot(ZC_range,PDF_ZC_all,'k')
hold