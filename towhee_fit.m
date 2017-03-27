function [ A, B, SIGA, SIGB ] = towhee_fit( X, Y, SIG )
%This is the fit subroutine in the towhee regression

WT = 1./(SIG.^2);
SS = sum(WT);
SX = sum(X.*WT);
SY = sum(Y.*WT);

SXOSS = SX/SS;

T = (X-SXOSS)./SIG;
ST2 = sum(T.*T);
B = sum(T.*Y./SIG);

B = B/ST2;
A = (SY - SX*B)/SS;
SIGA = sqrt((1+SX*SX/(SS*ST2))/SS);
SIGB = sqrt(1/ST2);

end


