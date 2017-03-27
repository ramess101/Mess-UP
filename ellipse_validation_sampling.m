% In this code I will sample:
% 20 points between 0 and 65
% 10 points between 65 and 85
% 50 points between 85 and 95
% 20 points between 95 and 99

% Must provide the PDF as a matrix titled "PDF" and epsilon and sigma as
% "epsilon" and "sigma"

clearvars -except PDF epsilon sigma

samples = [20; 10; 30; 40];

acceptable_epsilon = zeros(4,length(epsilon)*length(sigma));
acceptable_sigma = acceptable_epsilon;

k = ones(4,1);

for h = 1:length(epsilon)
    
    for i = 1:length(sigma)
        
        if PDF(h,i) <= 0.65
            
            acceptable_epsilon(1,k(1)) = epsilon(h);
            acceptable_sigma(1,k(1)) = sigma(i);
            
            k(1) = k(1) + 1;
            
        elseif PDF(h,i) <= 0.85
            
            acceptable_epsilon(2,k(2)) = epsilon(h);
            acceptable_sigma(2,k(2)) = sigma(i);
            
            k(2) = k(2) + 1;
            
        elseif PDF(h,i) <= 0.95

            acceptable_epsilon(3,k(3)) = epsilon(h);
            acceptable_sigma(3,k(3)) = sigma(i);
            
            k(3) = k(3) + 1;
            
        elseif PDF(h,i) <= 0.99
            
            acceptable_epsilon(4,k(4)) = epsilon(h);
            acceptable_sigma(4,k(4)) = sigma(i);
            
            k(4) = k(4) + 1;
                        
        end
        
    end
       
end

acceptable_65 = [acceptable_epsilon(1,:)]'; 
acceptable_85 = [acceptable_epsilon(2,:)]'; 
acceptable_95 = [acceptable_epsilon(3,:)]'; 
acceptable_99 = [acceptable_epsilon(4,:)]';

acceptable_65 = acceptable_65(acceptable_65~=0);
acceptable_85 = acceptable_85(acceptable_85~=0);
acceptable_95 = acceptable_95(acceptable_95~=0);
acceptable_99 = acceptable_99(acceptable_99~=0);

points = [length(acceptable_65),length(acceptable_85),length(acceptable_95),length(acceptable_99)];

acceptable_65(:,2) = [acceptable_sigma(1,1:points(1))]';
acceptable_85(:,2) = [acceptable_sigma(2,1:points(2))]';
acceptable_95(:,2) = [acceptable_sigma(3,1:points(3))]';
acceptable_99(:,2) = [acceptable_sigma(4,1:points(4))]';

figure
hold
scatter(acceptable_65(:,2),acceptable_65(:,1),'bx')
scatter(acceptable_85(:,2),acceptable_85(:,1),'rx')
scatter(acceptable_95(:,2),acceptable_95(:,1),'gx')
scatter(acceptable_99(:,2),acceptable_99(:,1),'cx')
hold

% Now we sample the amount prescribed from each region

k = 1;

for i = 1:4
    
    start(i) = k;
    
    for n = 1:samples(i)
        
        r = rand()*points(i);
        
        r = round(r);
        
        if r == 0
            
            r = 1;
            
        end
        
        eps(k) = acceptable_epsilon(i,r) + (1-rand())*0.01;
        sig(k) = acceptable_sigma(i,r) + (1-rand())*0.0001;
        
        k  = k + 1;
        
    end
    
    finish(i) = k - 1;
    
end

% figure
% hold
% scatter(sig(start(1):finish(1)),eps(start(1):finish(1)),'b')
% scatter(sig(start(2):finish(2)),eps(start(2):finish(2)),'r')
% scatter(sig(start(3):finish(3)),eps(start(3):finish(3)),'g')
% scatter(sig(start(4):finish(4)),eps(start(4):finish(4)),'k')
% xlim([3.966 3.983])
% ylim([45.20 45.55])
% hold

% hold
% scatter(sig,eps,'k')
% xlim([3.966 3.983])
% ylim([45.20 45.55])
% hold

% Without high temperature LDN
hold
scatter(sig,eps,'k')
xlim([3.963 3.981])
ylim([45.28 45.48])
hold
