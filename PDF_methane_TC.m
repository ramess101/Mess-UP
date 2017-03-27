clear

dTC = 0.00001;

TC_exp = 190.564;
TC_B2 = 195.0865;
TC_rhoL = 192.1065;

TC = [TC_B2, TC_rhoL, TC_exp];

n_TC = [0.120569, 0.119456];

p_TC = [0.2641, 0.035356];

sigma_exp = TC_exp * 0.002/3.29;

sigma_TC = sqrt(n_TC.^2+p_TC.^2);

sigma_TC(3) = sigma_exp;

low_TC = TC - 4*1.96*sigma_TC;

high_TC = TC + 4*1.96*sigma_TC;

for i = 1:3
    
    n = (high_TC - low_TC)./dTC;
    
    TC_i(i,1) = low_TC(i);
        
    for j = 2:n(i)
   
        TC_i(i,j) = TC_i(i,j-1) + dTC;
        
    end
    
end

TC_B2 = TC_i(1,:);
TC_B2 = TC_B2(TC_B2~=0);
TC_rhoL = TC_i(2,:);
TC_rhoL = TC_rhoL(TC_rhoL~=0);
TC_exp = TC_i(3,:);
TC_exp = TC_exp(TC_exp~=0);

n_B2 = pdf('Normal',TC_B2,TC(1),n_TC(1));
p_B2 = pdf('Normal',TC_B2,TC(1),p_TC(1));
overall_B2 = pdf('Normal',TC_B2,TC(1),sigma_TC(1));

n_rhoL = pdf('Normal',TC_rhoL,TC(2),n_TC(2));
p_rhoL = pdf('Normal',TC_rhoL,TC(2),p_TC(2));
overall_rhoL = pdf('Normal',TC_rhoL,TC(2),sigma_TC(2));

overall_exp = pdf('Normal',TC_exp,TC(3),sigma_TC(3));

n_B2 = n_B2/sum(n_B2);
p_B2 = p_B2/sum(p_B2);
overall_B2 = overall_B2/sum(overall_B2);

n_rhoL = n_rhoL/sum(n_rhoL);
p_rhoL = p_rhoL/sum(p_rhoL);
overall_rhoL = overall_rhoL/sum(overall_rhoL);

overall_exp = overall_exp/sum(overall_exp);

figure
hold
plot(TC_B2,n_B2,'r')
plot(TC_B2,p_B2,'b')
plot(TC_B2,overall_B2,'g')
hold

figure
hold
plot(TC_rhoL,n_rhoL,'r')
plot(TC_rhoL,p_rhoL,'b')
plot(TC_rhoL,overall_rhoL,'g')
hold

figure
plot(TC_exp,overall_exp,'k')

figure
hold
plot(TC_exp,overall_exp,'k-')
plot(TC_B2,n_B2,'r:')
plot(TC_B2,p_B2,'b:')
plot(TC_B2,overall_B2,'g:')
plot(TC_rhoL,n_rhoL,'r--')
plot(TC_rhoL,p_rhoL,'b--')
plot(TC_rhoL,overall_rhoL,'g--')
legend('Experimental','Numerical, B2','Parameter B2','Overall B2','Numerical rhoL','Parameter rhoL','Overall rhoL')
xlim([190 196])
xlabel('Critical Temperature (K)')
ylabel('PDF')
hold
