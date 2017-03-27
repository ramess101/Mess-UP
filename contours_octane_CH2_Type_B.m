clear

PDF_all=dlmread('PDF_octane_CH2_Type_B.txt');

eps_all = [45	45.0075758	45.0151515	45.0227273	45.0303030	45.0378788	45.0454545	45.0530303	45.0606061	45.0681818	45.0757576	45.0833333	45.0909091	45.0984848	45.1060606	45.1136364	45.1212121	45.1287879	45.1363636	45.1439394	45.1515152	45.1590909	45.1666667	45.1742424	45.1818182	45.1893939	45.1969697	45.2045455	45.2121212	45.2196970	45.2272727	45.2348485	45.2424242	45.2500000	45.2575758	45.2651515	45.2727273	45.2803030	45.2878788	45.2954545	45.3030303	45.3106061	45.3181818	45.3257576	45.3333333	45.3409091	45.3484848	45.3560606	45.3636364	45.3712121	45.3787879	45.3863636	45.3939394	45.4015152	45.4090909	45.4166667	45.4242424	45.4318182	45.4393939	45.4469697	45.4545455	45.4621212	45.4696970	45.4772727	45.4848485	45.4924242	45.5000000	45.5075758	45.5151515	45.5227273	45.5303030	45.5378788	45.5454545	45.5530303	45.5606061	45.5681818	45.5757576	45.5833333	45.5909091	45.5984848	45.6060606	45.6136364	45.6212121	45.6287879	45.6363636	45.6439394	45.6515152	45.6590909	45.6666667	45.6742424	45.6818182	45.6893939	45.6969697	45.7045455	45.7121212	45.7196970	45.7272727	45.7348485	45.7424242	45.7500000];
sig_all = [3.94	3.940707071	3.941414141	3.942121212	3.942828283	3.943535354	3.944242424	3.944949495	3.945656566	3.946363636	3.947070707	3.947777778	3.948484848	3.949191919	3.94989899	3.950606061	3.951313131	3.952020202	3.952727273	3.953434343	3.954141414	3.954848485	3.955555556	3.956262626	3.956969697	3.957676768	3.958383838	3.959090909	3.95979798	3.960505051	3.961212121	3.961919192	3.962626263	3.963333333	3.964040404	3.964747475	3.965454545	3.966161616	3.966868687	3.967575758	3.968282828	3.968989899	3.96969697	3.97040404	3.971111111	3.971818182	3.972525253	3.973232323	3.973939394	3.974646465	3.975353535	3.976060606	3.976767677	3.977474747	3.978181818	3.978888889	3.97959596	3.98030303	3.981010101	3.981717172	3.982424242	3.983131313	3.983838384	3.984545455	3.985252525	3.985959596	3.986666667	3.987373737	3.988080808	3.988787879	3.989494949	3.99020202	3.990909091	3.991616162	3.992323232	3.993030303	3.993737374	3.994444444	3.995151515	3.995858586	3.996565657	3.997272727	3.997979798	3.998686869	3.999393939	4.00010101	4.000808081	4.001515152	4.002222222	4.002929293	4.003636364	4.004343434	4.005050505	4.005757576	4.006464646	4.007171717	4.007878788	4.008585859	4.009292929	4.01];

s=1;

for i = 1:length(eps_all)
    
    for j = 1:length(sig_all)
    
        PDF(s) = PDF_all(i,j);
        eps(s) = eps_all(i);
        sig(s) = sig_all(j);

        s = s+1;

    end

end

eps_min = min(eps(PDF<=0.96));
eps_max = max(eps(PDF<=0.96));

sig_min = min(sig(PDF<=0.96));
sig_max = max(sig(PDF<=0.96));