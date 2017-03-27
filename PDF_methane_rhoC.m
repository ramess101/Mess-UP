clear

drhoC = 0.000002;

rhoC_exp = 10.142;
rhoC_B2 = 9.441018;
rhoC_rhoL = 10.23954;

rhoC = [rhoC_B2, rhoC_rhoL, rhoC_exp];

n_rhoC = [0.06003816, 0.06613];

p_rhoC = [0.04153878, 0.004764];

sigma_exp = rhoC_exp * 0.002/3.29;

sigma_rhoC = sqrt(n_rhoC.^2+p_rhoC.^2);

sigma_rhoC(3) = sigma_exp;

low_rhoC = rhoC - 4*1.96*sigma_rhoC;

high_rhoC = rhoC + 4*1.96*sigma_rhoC;

for i = 1:3
    
    n = (high_rhoC - low_rhoC)./drhoC;
    
    rhoC_i(i,1) = low_rhoC(i);
        
    for j = 2:n(i)
   
        rhoC_i(i,j) = rhoC_i(i,j-1) + drhoC;
        
    end
    
end

rhoC_B2 = rhoC_i(1,:);
rhoC_B2 = rhoC_B2(rhoC_B2~=0);
rhoC_rhoL = rhoC_i(2,:);
rhoC_rhoL = rhoC_rhoL(rhoC_rhoL~=0);
rhoC_exp = rhoC_i(3,:);
rhoC_exp = rhoC_exp(rhoC_exp~=0);

n_B2 = pdf('Normal',rhoC_B2,rhoC(1),n_rhoC(1));
p_B2 = pdf('Normal',rhoC_B2,rhoC(1),p_rhoC(1));
overall_B2 = pdf('Normal',rhoC_B2,rhoC(1),sigma_rhoC(1));

n_rhoL = pdf('Normal',rhoC_rhoL,rhoC(2),n_rhoC(2));
p_rhoL = pdf('Normal',rhoC_rhoL,rhoC(2),p_rhoC(2));
overall_rhoL = pdf('Normal',rhoC_rhoL,rhoC(2),sigma_rhoC(2));

overall_exp = pdf('Normal',rhoC_exp,rhoC(3),sigma_rhoC(3));

n_B2 = n_B2/sum(n_B2);
p_B2 = p_B2/sum(p_B2);
overall_B2 = overall_B2/sum(overall_B2);

n_rhoL = n_rhoL/sum(n_rhoL);
p_rhoL = p_rhoL/sum(p_rhoL);
overall_rhoL = overall_rhoL/sum(overall_rhoL);

overall_exp = overall_exp/sum(overall_exp);

figure
hold
plot(rhoC_B2,n_B2,'r')
plot(rhoC_B2,p_B2,'b')
plot(rhoC_B2,overall_B2,'g')
hold

figure
hold
plot(rhoC_rhoL,n_rhoL,'r')
plot(rhoC_rhoL,p_rhoL,'b')
plot(rhoC_rhoL,overall_rhoL,'g')
hold

figure
plot(rhoC_exp,overall_exp,'k')

figure
hold
plot(rhoC_exp,overall_exp,'k-')
plot(rhoC_B2,n_B2,'r:')
plot(rhoC_B2,p_B2,'b:')
plot(rhoC_B2,overall_B2,'g:')
plot(rhoC_rhoL,n_rhoL,'r--')
plot(rhoC_rhoL,p_rhoL,'b--')
plot(rhoC_rhoL,overall_rhoL,'g--')
legend('Experimental','Numerical, B2','Parameter B2','Overall B2','Numerical rhoL','Parameter rhoL','Overall rhoL')
xlim([9.2 10.4])
xlabel('Critical Density (mol/L)')
ylabel('PDF')
hold
