% PDF_all = zeros(length(epsilon_refined)*length(sigma_refined),6+length(Temp_exp));

PDF_all = zeros(length(epsilon_refined)*length(sigma_refined),6);
index = 1;

for i = 1:length(epsilon_refined)
    
    for j = 1:length(sigma_refined)
        
        PDF_all(index,1) = PDF(i,j);
        PDF_all(index,2) = epsilon_refined(i);
        PDF_all(index,3) = sigma_refined(j);
        PDF_all(index,4) = TC_opt(i,j);
        PDF_all(index,5) = rhoc_opt(i,j);
        PDF_all(index,6) = PC(i,j);
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
rhoc_occ = PDF_occ;
PC_occ = PDF_occ;
% rhoL_occ = zeros(length(Temp_exp),length(PDF_centers),max(PDF_counts));

for i = 1:sum(PDF_counts)
       
    dif_centers = abs(PDF_all(i,1) - PDF_centers);
    [~,I] = min(dif_centers);
    PDF_counted(I) = PDF_counted(I) + 1;
    PDF_occ(I,PDF_counted(I)) = PDF_all(i,1);
    eps_occ(I,PDF_counted(I)) = PDF_all(i,2);
    sig_occ(I,PDF_counted(I)) = PDF_all(i,3);
    TC_occ(I,PDF_counted(I)) = PDF_all(i,4);
    rhoc_occ(I,PDF_counted(I)) = PDF_all(i,5);
    PC_occ(I,PDF_counted(I)) = PDF_all(i,6);
    
%         for k = 1:length(Temp_exp)
%     
%             rhoL_occ(k,I,PDF_counted(I)) = PDF_all(i,6+k);
%     
%         end
    
end

ii = 1000000;

eps_norm = zeros(ii,1);
sig_norm = eps_norm;
TC_norm = eps_norm;
rhoc_norm = eps_norm;
PC_norm = eps_norm;
% rhoL_norm = zeros(ii,length(Temp_exp));

for i = 1:ii
   
    r=fcdf(frnd(p,n-p),p,n-p);
    
    % I think a better way is to organize the PDF_centers and bin them such
    % that there are ~10 eps/sig in each bin. Then this new discrete one
    % dimensional PDF can be sampled from and the choice of which of the 10
    % eps/sig is assigned would be random. Need to keep track of the
    % indices. Also probably need to convert NxN grid to N^2 column
    
    dif = abs(r-PDF_centers);
    
    [dif_min, I] = min(dif);
    
    J = ceil(rand()*PDF_counts(I));
    
    eps_norm(i) = eps_occ(I,J);
    
    sig_norm(i) = sig_occ(I,J);
    
    TC_norm(i) = TC_occ(I,J);
    
    rhoc_norm(i) = rhoc_occ(I,J);
    
    PC_norm(i) = PC_occ(I,J);
%     
%         for k = 1:length(Temp_exp)
%     
%             rhoL_norm(i,k) = rhoL_occ(k,I,J);
%     
%         end
end

ZC_norm = PC_norm * MW ./ rhoc_norm ./ TC_norm / 8.314472;

clear eps_occ sig_occ TC_occ rhoc_occ PC_occ rhoL_occ
