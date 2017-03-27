clear

dPC = 0.000001;

PC_exp = 4.599;
PC_B2 = 4.702777;
PC_rhoL = 5.022624;

PC = [PC_B2, PC_rhoL, PC_exp];

n_PC = [0.02448343, 0.026149];

p_PC = [0.02780562, 0.001853];

sigma_exp = PC_exp * 0.002/3.29;

sigma_PC = sqrt(n_PC.^2+p_PC.^2);

sigma_PC(3) = sigma_exp;

low_PC = PC - 4*1.96*sigma_PC;

high_PC = PC + 4*1.96*sigma_PC;

for i = 1:3
    
    n = (high_PC - low_PC)./dPC;
    
    PC_i(i,1) = low_PC(i);
        
    for j = 2:n(i)
   
        PC_i(i,j) = PC_i(i,j-1) + dPC;
        
    end
    
end

PC_B2 = PC_i(1,:);
PC_B2 = PC_B2(PC_B2~=0);
PC_rhoL = PC_i(2,:);
PC_rhoL = PC_rhoL(PC_rhoL~=0);
PC_exp = PC_i(3,:);
PC_exp = PC_exp(PC_exp~=0);

n_B2 = pdf('Normal',PC_B2,PC(1),n_PC(1));
p_B2 = pdf('Normal',PC_B2,PC(1),p_PC(1));
overall_B2 = pdf('Normal',PC_B2,PC(1),sigma_PC(1));

n_rhoL = pdf('Normal',PC_rhoL,PC(2),n_PC(2));
p_rhoL = pdf('Normal',PC_rhoL,PC(2),p_PC(2));
overall_rhoL = pdf('Normal',PC_rhoL,PC(2),sigma_PC(2));

overall_exp = pdf('Normal',PC_exp,PC(3),sigma_PC(3));

n_B2 = n_B2/sum(n_B2);
p_B2 = p_B2/sum(p_B2);
overall_B2 = overall_B2/sum(overall_B2);

n_rhoL = n_rhoL/sum(n_rhoL);
p_rhoL = p_rhoL/sum(p_rhoL);
overall_rhoL = overall_rhoL/sum(overall_rhoL);

overall_exp = overall_exp/sum(overall_exp);

figure
hold
plot(PC_B2,n_B2,'r')
plot(PC_B2,p_B2,'b')
plot(PC_B2,overall_B2,'g')
hold

figure
hold
plot(PC_rhoL,n_rhoL,'r')
plot(PC_rhoL,p_rhoL,'b')
plot(PC_rhoL,overall_rhoL,'g')
hold

figure
plot(PC_exp,overall_exp,'k')

figure
hold
plot(PC_exp,overall_exp,'k-')
plot(PC_B2,n_B2,'r:')
plot(PC_B2,p_B2,'b:')
plot(PC_B2,overall_B2,'g:')
plot(PC_rhoL,n_rhoL,'r--')
plot(PC_rhoL,p_rhoL,'b--')
plot(PC_rhoL,overall_rhoL,'g--')
legend('Experimental','Numerical, B2','Parameter B2','Overall B2','Numerical rhoL','Parameter rhoL','Overall rhoL')
xlim([4.5 5.1])
xlabel('Critical Pressure (MPa)')
ylabel('PDF')
hold
