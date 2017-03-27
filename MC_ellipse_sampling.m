PDF_all = zeros(length(epsilon)*length(sigma),3);
index = 1;

p = 2;
n = 11; % 12 for ethane, 11 for octane

for i = 1:length(epsilon)
    
    for j = 1:length(sigma)
        
        PDF_all(index,1) = PDF(i,j);
        PDF_all(index,2) = epsilon(i);
        PDF_all(index,3) = sigma(j);
                
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

for i = 1:sum(PDF_counts)
   
    dif_centers = abs(PDF_all(i,1) - PDF_centers);
    [~,I] = min(dif_centers);
    PDF_counted(I) = PDF_counted(I) + 1;
    PDF_occ(I,PDF_counted(I)) = PDF_all(i,1);
    eps_occ(I,PDF_counted(I)) = PDF_all(i,2);
    sig_occ(I,PDF_counted(I)) = PDF_all(i,3);
        
end

ii = 20000;

eps_norm = zeros(ii,1);
sig_norm = eps_norm;

for i = 1:ii
   
    r=fcdf(frnd(p,n-p),p,n-p);
    
    dif = abs(r-PDF_centers);
    
    [dif_min, I] = min(dif);
    
    J = ceil(rand()*PDF_counts(I));
    
    eps_norm(i) = eps_occ(I,J);
    
    sig_norm(i) = sig_occ(I,J);
    
end

clear eps_occ sig_occ TC_occ rhoc_occ PC_occ rhoL_occ