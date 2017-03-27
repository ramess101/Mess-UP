% Need to run the MC file first to get all the vap_ave and liq_ave

% TC_range = linspace(303.75,305.25,30); % For Type A

TC_range = linspace(300,308.5,170); % For Type AB

PDF_TC = zeros(length(sigma),length(TC_range));

for m = 56:length(sigma)
    
    T = Temp;
    pv = vap_ave(:,m)';
    pL = liq_ave(:,m)';
       
[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma2 = (SSE_fit/(n-p))^2;

A_range = linspace(0.9775,1.0225,30)*A_fit;
b_range = linspace(0.9975,1.0025,30)*b_fit;
pc_range = linspace(0.995,1.005,30)*pc_fit;

for g=1:length(A_range)
    
    for h=1:length(b_range)
    
        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
        
        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
        
        SSE = sum(SE);
        
        PDF_TC(m,j) = PDF_TC(m,j) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);

            end 
            
        end
    
    end
    
end

PDF_TC(m,:) = PDF_TC(m,:)/sum(PDF_TC(m,:));

end

PDF_TC_all = sum(PDF_TC);

PDF_TC_all = PDF_TC_all/sum(PDF_TC_all);

figure
hold
plot(TC_range,PDF_TC_all,'k')
hold