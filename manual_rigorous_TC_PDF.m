[TC_fit, pc_fit, A_fit, b_fit, sy, sz] = towhee_error_model(T,pv,pL,beta);

n = 2*length(T);

p = 4;

y = (pv + pL)/2;
z = (pL - pv)/2;

SE_fit = ((y - (pc_fit + A_fit*(TC_fit-T)))./sy).^2 + ((z - (b_fit*(TC_fit-T).^beta))./sz).^2; % This includes uncertainties

SSE_fit = sum(SE_fit);

sigma2 = (SSE_fit/(n-p))^2;

% Type A analysis

% A_range = linspace(0.9775,1.0225,20)*A_fit;
% b_range = linspace(0.9975,1.0025,20)*b_fit;
% pc_range = linspace(0.995,1.005,20)*pc_fit;
% TC_range = linspace(0.997,1.003,20)*TC_fit;

% For MCO

% A_range = linspace(0.93,1.07,20)*A_fit;
% b_range = linspace(0.99,1.01,20)*b_fit;
% pc_range = linspace(0.98,1.02,20)*pc_fit;
% TC_range = linspace(0.995,1.005,20)*TC_fit;

% For MC_C16

% A_range = linspace(0.96,1.04,20)*A_fit;
% b_range = linspace(0.995,1.005,20)*b_fit;
% pc_range = linspace(0.985,1.015,20)*pc_fit;
% TC_range = linspace(0.995,1.005,20)*TC_fit;

% For second parameter set MC_C16
% A_range = linspace(0.96,1.04,30)*A_fit;
% b_range = linspace(0.995,1.005,30)*b_fit;
% pc_range = linspace(0.985,1.015,30)*pc_fit;
% TC_range = linspace(726,728,40);

% For MC_C16 with Type AB

A_range = linspace(0.96,1.04,30)*A_fit;
b_range = linspace(0.995,1.005,30)*b_fit;
pc_range = linspace(0.985,1.015,30)*pc_fit;
TC_range = linspace(727,733,40);

% For MC_C24

% A_range = linspace(0.9,1.1,20)*A_fit;
% b_range = linspace(0.98,1.02,20)*b_fit;
% pc_range = linspace(0.95,1.05,20)*pc_fit;
% TC_range = linspace(0.99,1.01,20)*TC_fit;

% For MC_C36

% A_range = linspace(0.75,1.25,20)*A_fit;
% b_range = linspace(0.96,1.04,20)*b_fit;
% pc_range = linspace(0.85,1.15,20)*pc_fit;
% TC_range = linspace(0.975,1.025,20)*TC_fit;

% For MC_C48

% A_range = linspace(0.75,1.25,20)*A_fit;
% b_range = linspace(0.96,1.04,20)*b_fit;
% pc_range = linspace(0.85,1.15,20)*pc_fit;
% TC_range = linspace(0.975,1.025,20)*TC_fit;

PDF_A = zeros(length(A_range),1);
PDF_b = zeros(length(b_range),1);
PDF_pc = zeros(length(pc_range),1);
PDF_TC = zeros(length(TC_range),1);

for g=1:length(A_range)
    
    for h=1:length(b_range)
    
        for i=1:length(pc_range)
    
            for j=1:length(TC_range)
        
        SE = ((y - (pc_range(i) + A_range(g)*(TC_range(j)-T)))./sy).^2 + ((z - (b_range(h)*(TC_range(j)-T).^beta))./sz).^2;
        
        SSE = sum(SE);
        
        PDF_A(g) = PDF_A(g) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
        
        PDF_b(h) = PDF_b(h) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
        
        PDF_pc(i) = PDF_pc(i) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p); 
        
        PDF_TC(j) = PDF_TC(j) + fpdf((SSE - SSE_fit)/(sigma2*p),p,n-p);
        

            end 
            
        end
    
    end
    
end

PDF_A = PDF_A / sum(PDF_A);
PDF_b = PDF_b / sum(PDF_b);
PDF_pc = PDF_pc / sum(PDF_pc);
PDF_TC = PDF_TC / sum(PDF_TC);

figure
subplot(2,2,1)
scatter(A_range,PDF_A)
subplot(2,2,2)
scatter(b_range,PDF_b)
subplot(2,2,3)
scatter(pc_range,PDF_pc)
subplot(2,2,4)
scatter(TC_range,PDF_TC)

