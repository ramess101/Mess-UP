clear

LMOL = [21.61 21.32 20.45 19.47 18.52 17.42 16.3 14.89 13.24 10.1 8.82]; %mol/L

Temp = [92 100 124 150 174 200 224 250 274 300 304];

MW = 30.06904; %gm/mol

LDN = LMOL*MW/1000; % gm/mL

Rack = [1.9122, 0.27937, 305.32, 0.29187];

LMOL_Rack = Rack(1)./(Rack(2).^(1+(1-Temp/Rack(3)).^Rack(4))); %mol/L

LDN_Rack = LMOL_Rack * MW / 1000; %gm/mL

SSE = sum((LDN-LDN_Rack).^2);

n = length(Temp);

p = length(Rack);

sigma = sqrt(SSE/(n-p));