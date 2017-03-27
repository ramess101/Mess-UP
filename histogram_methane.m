% PDF_all = zeros(length(epsilon_refined)*length(sigma_refined),6+length(Temp_exp));

PDF_all = zeros(length(eps_range)*length(sig_range),4);
index = 1;

for i = 1:length(eps_range)
    
    for j = 1:length(sig_range)
        
        PDF_all(index,1) = PDF(i,j);
        PDF_all(index,2) = eps_range(i);
        PDF_all(index,3) = sig_range(j);
        PDF_all(index,4) = TC_range(i,j);
%         PDF_all(index,5) = rhoC_range(i,j);
%         PDF_all(index,6) = PC_range(i,j);
        
%         PDF_all(index,2) = eps_sig3_range(i,j);
%         
%             for k = 1:length(Temp_exp)
%                 
%              PDF_all(index,6+k) = liq_hat(k,i,j);
%         
%             end
        
        index = index + 1;
        
    end
    
end

[PDF_counts, PDF_centers] = hist(PDF_all(:,1),500);

PDF_centers = PDF_centers(PDF_counts~=0);
PDF_counts = PDF_counts(PDF_counts~=0);
PDF_counted = zeros(length(PDF_counts),1);
PDF_occ = zeros(length(PDF_centers),max(PDF_counts));
eps_occ = PDF_occ;
sig_occ = PDF_occ;
TC_occ = PDF_occ;
% rhoC_occ = PDF_occ;
% PC_occ = PDF_occ;
% eps_sig3_occ = PDF_occ;
% rhoL_occ = zeros(length(Temp_exp),length(PDF_centers),max(PDF_counts));

for i = 1:sum(PDF_counts)
   
    dif_centers = abs(PDF_all(i,1) - PDF_centers);
    [~,I] = min(dif_centers);
    PDF_counted(I) = PDF_counted(I) + 1;
    PDF_occ(I,PDF_counted(I)) = PDF_all(i,1);
    eps_occ(I,PDF_counted(I)) = PDF_all(i,2);
    sig_occ(I,PDF_counted(I)) = PDF_all(i,3);
    TC_occ(I,PDF_counted(I)) = PDF_all(i,4);
%     rhoC_occ(I,PDF_counted(I)) = PDF_all(i,5);
%     PC_occ(I,PDF_counted(I)) = PDF_all(i,6);
%     eps_sig3_occ(I,PDF_counted(I)) = PDF_all(i,2);
%     
%         for k = 1:length(Temp_exp)
%     
%             rhoL_occ(k,I,PDF_counted(I)) = PDF_all(i,6+k);
%     
%         end
    
end

ii = 500000;

eps_norm = zeros(ii,1);
sig_norm = eps_norm;
TC_norm = eps_norm;
% rhoC_norm = eps_norm;
% PC_norm = eps_norm;
% eps_sig3_norm = eps_norm;
% rhoL_norm = zeros(ii,length(Temp_exp));

for i = 1:ii
   
    r=fcdf(frnd(p,n-p),p,n-p);
    
    dif = abs(r-PDF_centers);
    
    [dif_min, I] = min(dif);
    
    J = ceil(rand()*PDF_counts(I));
    
    eps_norm(i) = eps_occ(I,J);
    
    sig_norm(i) = sig_occ(I,J);
    
    TC_norm(i) = TC_occ(I,J);
%     
%     rhoC_norm(i) = rhoC_occ(I,J);
%     
%     PC_norm(i) = PC_occ(I,J);

%     eps_sig3_norm(i) = eps_sig3_occ(I,J);
%     
%         for k = 1:length(Temp_exp)
%     
%             rhoL_norm(i,k) = rhoL_occ(k,I,J);
%     
%         end
end

clear eps_occ sig_occ TC_occ rhoC_occ PC_occ rhoL_occ
