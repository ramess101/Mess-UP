clear

PDF_all=dlmread('PDF_methane_rhoL_coarse.txt');

% 200x200
% eps_all = [146.1870459	146.1885166	146.1899873	146.191458	146.1929287	146.1943994	146.1958701	146.1973407	146.1988114	146.2002821	146.2017528	146.2032235	146.2046942	146.2061649	146.2076356	146.2091062	146.2105769	146.2120476	146.2135183	146.214989	146.2164597	146.2179304	146.2194011	146.2208717	146.2223424	146.2238131	146.2252838	146.2267545	146.2282252	146.2296959	146.2311665	146.2326372	146.2341079	146.2355786	146.2370493	146.23852	146.2399907	146.2414614	146.242932	146.2444027	146.2458734	146.2473441	146.2488148	146.2502855	146.2517562	146.2532269	146.2546975	146.2561682	146.2576389	146.2591096	146.2605803	146.262051	146.2635217	146.2649924	146.266463	146.2679337	146.2694044	146.2708751	146.2723458	146.2738165	146.2752872	146.2767579	146.2782285	146.2796992	146.2811699	146.2826406	146.2841113	146.285582	146.2870527	146.2885234	146.289994	146.2914647	146.2929354	146.2944061	146.2958768	146.2973475	146.2988182	146.3002888	146.3017595	146.3032302	146.3047009	146.3061716	146.3076423	146.309113	146.3105837	146.3120543	146.313525	146.3149957	146.3164664	146.3179371	146.3194078	146.3208785	146.3223492	146.3238198	146.3252905	146.3267612	146.3282319	146.3297026	146.3311733	146.332644	146.3341147	146.3355853	146.337056	146.3385267	146.3399974	146.3414681	146.3429388	146.3444095	146.3458802	146.3473508	146.3488215	146.3502922	146.3517629	146.3532336	146.3547043	146.356175	146.3576456	146.3591163	146.360587	146.3620577	146.3635284	146.3649991	146.3664698	146.3679405	146.3694111	146.3708818	146.3723525	146.3738232	146.3752939	146.3767646	146.3782353	146.379706	146.3811766	146.3826473	146.384118	146.3855887	146.3870594	146.3885301	146.3900008	146.3914715	146.3929421	146.3944128	146.3958835	146.3973542	146.3988249	146.4002956	146.4017663	146.403237	146.4047076	146.4061783	146.407649	146.4091197	146.4105904	146.4120611	146.4135318	146.4150025	146.4164731	146.4179438	146.4194145	146.4208852	146.4223559	146.4238266	146.4252973	146.4267679	146.4282386	146.4297093	146.43118	146.4326507	146.4341214	146.4355921	146.4370628	146.4385334	146.4400041	146.4414748	146.4429455	146.4444162	146.4458869	146.4473576	146.4488283	146.4502989	146.4517696	146.4532403	146.454711	146.4561817	146.4576524	146.4591231	146.4605938	146.4620644	146.4635351	146.4650058	146.4664765	146.4679472	146.4694179	146.4708886	146.4723593	146.4738299	146.4753006	146.4767713	146.478242	146.4797127];
% sig_all = [3.710686073	3.710723403	3.710760734	3.710798065	3.710835395	3.710872726	3.710910057	3.710947387	3.710984718	3.711022049	3.711059379	3.71109671	3.71113404	3.711171371	3.711208702	3.711246032	3.711283363	3.711320694	3.711358024	3.711395355	3.711432686	3.711470016	3.711507347	3.711544678	3.711582008	3.711619339	3.71165667	3.711694	3.711731331	3.711768662	3.711805992	3.711843323	3.711880654	3.711917984	3.711955315	3.711992646	3.712029976	3.712067307	3.712104638	3.712141968	3.712179299	3.71221663	3.71225396	3.712291291	3.712328622	3.712365952	3.712403283	3.712440614	3.712477944	3.712515275	3.712552606	3.712589936	3.712627267	3.712664597	3.712701928	3.712739259	3.712776589	3.71281392	3.712851251	3.712888581	3.712925912	3.712963243	3.713000573	3.713037904	3.713075235	3.713112565	3.713149896	3.713187227	3.713224557	3.713261888	3.713299219	3.713336549	3.71337388	3.713411211	3.713448541	3.713485872	3.713523203	3.713560533	3.713597864	3.713635195	3.713672525	3.713709856	3.713747187	3.713784517	3.713821848	3.713859179	3.713896509	3.71393384	3.713971171	3.714008501	3.714045832	3.714083162	3.714120493	3.714157824	3.714195154	3.714232485	3.714269816	3.714307146	3.714344477	3.714381808	3.714419138	3.714456469	3.7144938	3.71453113	3.714568461	3.714605792	3.714643122	3.714680453	3.714717784	3.714755114	3.714792445	3.714829776	3.714867106	3.714904437	3.714941768	3.714979098	3.715016429	3.71505376	3.71509109	3.715128421	3.715165752	3.715203082	3.715240413	3.715277744	3.715315074	3.715352405	3.715389736	3.715427066	3.715464397	3.715501727	3.715539058	3.715576389	3.715613719	3.71565105	3.715688381	3.715725711	3.715763042	3.715800373	3.715837703	3.715875034	3.715912365	3.715949695	3.715987026	3.716024357	3.716061687	3.716099018	3.716136349	3.716173679	3.71621101	3.716248341	3.716285671	3.716323002	3.716360333	3.716397663	3.716434994	3.716472325	3.716509655	3.716546986	3.716584317	3.716621647	3.716658978	3.716696309	3.716733639	3.71677097	3.716808301	3.716845631	3.716882962	3.716920292	3.716957623	3.716994954	3.717032284	3.717069615	3.717106946	3.717144276	3.717181607	3.717218938	3.717256268	3.717293599	3.71733093	3.71736826	3.717405591	3.717442922	3.717480252	3.717517583	3.717554914	3.717592244	3.717629575	3.717666906	3.717704236	3.717741567	3.717778898	3.717816228	3.717853559	3.71789089	3.71792822	3.717965551	3.718002882	3.718040212	3.718077543	3.718114874];

% 500x500
% eps_all = [146.187045931271,146.187632437801,146.188218944332,146.188805450862,146.189391957392,146.189978463922,146.190564970453,146.191151476983,146.191737983513,146.192324490044,146.192910996574,146.193497503104,146.194084009635,146.194670516165,146.195257022695,146.195843529225,146.196430035756,146.197016542286,146.197603048816,146.198189555347,146.198776061877,146.199362568407,146.199949074938,146.200535581468,146.201122087998,146.201708594528,146.202295101059,146.202881607589,146.203468114119,146.204054620650,146.204641127180,146.205227633710,146.205814140241,146.206400646771,146.206987153301,146.207573659832,146.208160166362,146.208746672892,146.209333179422,146.209919685953,146.210506192483,146.211092699013,146.211679205544,146.212265712074,146.212852218604,146.213438725135,146.214025231665,146.214611738195,146.215198244725,146.215784751256,146.216371257786,146.216957764316,146.217544270847,146.218130777377,146.218717283907,146.219303790438,146.219890296968,146.220476803498,146.221063310028,146.221649816559,146.222236323089,146.222822829619,146.223409336150,146.223995842680,146.224582349210,146.225168855741,146.225755362271,146.226341868801,146.226928375332,146.227514881862,146.228101388392,146.228687894922,146.229274401453,146.229860907983,146.230447414513,146.231033921044,146.231620427574,146.232206934104,146.232793440635,146.233379947165,146.233966453695,146.234552960225,146.235139466756,146.235725973286,146.236312479816,146.236898986347,146.237485492877,146.238071999407,146.238658505938,146.239245012468,146.239831518998,146.240418025528,146.241004532059,146.241591038589,146.242177545119,146.242764051650,146.243350558180,146.243937064710,146.244523571241,146.245110077771,146.245696584301,146.246283090832,146.246869597362,146.247456103892,146.248042610422,146.248629116953,146.249215623483,146.249802130013,146.250388636544,146.250975143074,146.251561649604,146.252148156135,146.252734662665,146.253321169195,146.253907675725,146.254494182256,146.255080688786,146.255667195316,146.256253701847,146.256840208377,146.257426714907,146.258013221438,146.258599727968,146.259186234498,146.259772741028,146.260359247559,146.260945754089,146.261532260619,146.262118767150,146.262705273680,146.263291780210,146.263878286741,146.264464793271,146.265051299801,146.265637806331,146.266224312862,146.266810819392,146.267397325922,146.267983832453,146.268570338983,146.269156845513,146.269743352044,146.270329858574,146.270916365104,146.271502871635,146.272089378165,146.272675884695,146.273262391225,146.273848897756,146.274435404286,146.275021910816,146.275608417347,146.276194923877,146.276781430407,146.277367936938,146.277954443468,146.278540949998,146.279127456528,146.279713963059,146.280300469589,146.280886976119,146.281473482650,146.282059989180,146.282646495710,146.283233002241,146.283819508771,146.284406015301,146.284992521831,146.285579028362,146.286165534892,146.286752041422,146.287338547953,146.287925054483,146.288511561013,146.289098067544,146.289684574074,146.290271080604,146.290857587135,146.291444093665,146.292030600195,146.292617106725,146.293203613256,146.293790119786,146.294376626316,146.294963132847,146.295549639377,146.296136145907,146.296722652438,146.297309158968,146.297895665498,146.298482172028,146.299068678559,146.299655185089,146.300241691619,146.300828198150,146.301414704680,146.302001211210,146.302587717741,146.303174224271,146.303760730801,146.304347237331,146.304933743862,146.305520250392,146.306106756922,146.306693263453,146.307279769983,146.307866276513,146.308452783044,146.309039289574,146.309625796104,146.310212302635,146.310798809165,146.311385315695,146.311971822225,146.312558328756,146.313144835286,146.313731341816,146.314317848347,146.314904354877,146.315490861407,146.316077367938,146.316663874468,146.317250380998,146.317836887528,146.318423394059,146.319009900589,146.319596407119,146.320182913650,146.320769420180,146.321355926710,146.321942433241,146.322528939771,146.323115446301,146.323701952831,146.324288459362,146.324874965892,146.325461472422,146.326047978953,146.326634485483,146.327220992013,146.327807498544,146.328394005074,146.328980511604,146.329567018135,146.330153524665,146.330740031195,146.331326537725,146.331913044256,146.332499550786,146.333086057316,146.333672563847,146.334259070377,146.334845576907,146.335432083438,146.336018589968,146.336605096498,146.337191603028,146.337778109559,146.338364616089,146.338951122619,146.339537629150,146.340124135680,146.340710642210,146.341297148741,146.341883655271,146.342470161801,146.343056668331,146.343643174862,146.344229681392,146.344816187922,146.345402694453,146.345989200983,146.346575707513,146.347162214044,146.347748720574,146.348335227104,146.348921733635,146.349508240165,146.350094746695,146.350681253225,146.351267759756,146.351854266286,146.352440772816,146.353027279347,146.353613785877,146.354200292407,146.354786798938,146.355373305468,146.355959811998,146.356546318528,146.357132825059,146.357719331589,146.358305838119,146.358892344650,146.359478851180,146.360065357710,146.360651864241,146.361238370771,146.361824877301,146.362411383831,146.362997890362,146.363584396892,146.364170903422,146.364757409953,146.365343916483,146.365930423013,146.366516929544,146.367103436074,146.367689942604,146.368276449134,146.368862955665,146.369449462195,146.370035968725,146.370622475256,146.371208981786,146.371795488316,146.372381994847,146.372968501377,146.373555007907,146.374141514438,146.374728020968,146.375314527498,146.375901034028,146.376487540559,146.377074047089,146.377660553619,146.378247060150,146.378833566680,146.379420073210,146.380006579741,146.380593086271,146.381179592801,146.381766099331,146.382352605862,146.382939112392,146.383525618922,146.384112125453,146.384698631983,146.385285138513,146.385871645044,146.386458151574,146.387044658104,146.387631164634,146.388217671165,146.388804177695,146.389390684225,146.389977190756,146.390563697286,146.391150203816,146.391736710347,146.392323216877,146.392909723407,146.393496229938,146.394082736468,146.394669242998,146.395255749528,146.395842256059,146.396428762589,146.397015269119,146.397601775650,146.398188282180,146.398774788710,146.399361295241,146.399947801771,146.400534308301,146.401120814831,146.401707321362,146.402293827892,146.402880334422,146.403466840953,146.404053347483,146.404639854013,146.405226360544,146.405812867074,146.406399373604,146.406985880134,146.407572386665,146.408158893195,146.408745399725,146.409331906256,146.409918412786,146.410504919316,146.411091425847,146.411677932377,146.412264438907,146.412850945438,146.413437451968,146.414023958498,146.414610465028,146.415196971559,146.415783478089,146.416369984619,146.416956491150,146.417542997680,146.418129504210,146.418716010741,146.419302517271,146.419889023801,146.420475530331,146.421062036862,146.421648543392,146.422235049922,146.422821556453,146.423408062983,146.423994569513,146.424581076044,146.425167582574,146.425754089104,146.426340595634,146.426927102165,146.427513608695,146.428100115225,146.428686621756,146.429273128286,146.429859634816,146.430446141347,146.431032647877,146.431619154407,146.432205660938,146.432792167468,146.433378673998,146.433965180528,146.434551687059,146.435138193589,146.435724700119,146.436311206650,146.436897713180,146.437484219710,146.438070726241,146.438657232771,146.439243739301,146.439830245831,146.440416752362,146.441003258892,146.441589765422,146.442176271953,146.442762778483,146.443349285013,146.443935791544,146.444522298074,146.445108804604,146.445695311134,146.446281817665,146.446868324195,146.447454830725,146.448041337256,146.448627843786,146.449214350316,146.449800856847,146.450387363377,146.450973869907,146.451560376438,146.452146882968,146.452733389498,146.453319896028,146.453906402559,146.454492909089,146.455079415619,146.455665922150,146.456252428680,146.456838935210,146.457425441741,146.458011948271,146.458598454801,146.459184961331,146.459771467862,146.460357974392,146.460944480922,146.461530987453,146.462117493983,146.462704000513,146.463290507044,146.463877013574,146.464463520104,146.465050026635,146.465636533165,146.466223039695,146.466809546225,146.467396052756,146.467982559286,146.468569065816,146.469155572347,146.469742078877,146.470328585407,146.470915091938,146.471501598468,146.472088104998,146.472674611528,146.473261118059,146.473847624589,146.474434131119,146.475020637650,146.475607144180,146.476193650710,146.476780157241,146.477366663771,146.477953170301,146.478539676831,146.479126183362,146.479712689892];
% sig_all = [3.71068607260122,3.71070095997787,3.71071584735451,3.71073073473116,3.71074562210780,3.71076050948445,3.71077539686110,3.71079028423774,3.71080517161439,3.71082005899103,3.71083494636768,3.71084983374432,3.71086472112097,3.71087960849761,3.71089449587426,3.71090938325091,3.71092427062755,3.71093915800420,3.71095404538084,3.71096893275749,3.71098382013413,3.71099870751078,3.71101359488742,3.71102848226407,3.71104336964072,3.71105825701736,3.71107314439401,3.71108803177065,3.71110291914730,3.71111780652394,3.71113269390059,3.71114758127724,3.71116246865388,3.71117735603053,3.71119224340717,3.71120713078382,3.71122201816046,3.71123690553711,3.71125179291375,3.71126668029040,3.71128156766705,3.71129645504369,3.71131134242034,3.71132622979698,3.71134111717363,3.71135600455027,3.71137089192692,3.71138577930356,3.71140066668021,3.71141555405686,3.71143044143350,3.71144532881015,3.71146021618679,3.71147510356344,3.71148999094008,3.71150487831673,3.71151976569338,3.71153465307002,3.71154954044667,3.71156442782331,3.71157931519996,3.71159420257660,3.71160908995325,3.71162397732989,3.71163886470654,3.71165375208319,3.71166863945983,3.71168352683648,3.71169841421312,3.71171330158977,3.71172818896641,3.71174307634306,3.71175796371970,3.71177285109635,3.71178773847300,3.71180262584964,3.71181751322629,3.71183240060293,3.71184728797958,3.71186217535622,3.71187706273287,3.71189195010951,3.71190683748616,3.71192172486281,3.71193661223945,3.71195149961610,3.71196638699274,3.71198127436939,3.71199616174603,3.71201104912268,3.71202593649933,3.71204082387597,3.71205571125262,3.71207059862926,3.71208548600591,3.71210037338255,3.71211526075920,3.71213014813584,3.71214503551249,3.71215992288914,3.71217481026578,3.71218969764243,3.71220458501907,3.71221947239572,3.71223435977236,3.71224924714901,3.71226413452565,3.71227902190230,3.71229390927895,3.71230879665559,3.71232368403224,3.71233857140888,3.71235345878553,3.71236834616217,3.71238323353882,3.71239812091546,3.71241300829211,3.71242789566876,3.71244278304540,3.71245767042205,3.71247255779869,3.71248744517534,3.71250233255198,3.71251721992863,3.71253210730528,3.71254699468192,3.71256188205857,3.71257676943521,3.71259165681186,3.71260654418850,3.71262143156515,3.71263631894179,3.71265120631844,3.71266609369508,3.71268098107173,3.71269586844838,3.71271075582502,3.71272564320167,3.71274053057831,3.71275541795496,3.71277030533160,3.71278519270825,3.71280008008490,3.71281496746154,3.71282985483819,3.71284474221483,3.71285962959148,3.71287451696812,3.71288940434477,3.71290429172141,3.71291917909806,3.71293406647471,3.71294895385135,3.71296384122800,3.71297872860464,3.71299361598129,3.71300850335793,3.71302339073458,3.71303827811122,3.71305316548787,3.71306805286452,3.71308294024116,3.71309782761781,3.71311271499445,3.71312760237110,3.71314248974774,3.71315737712439,3.71317226450104,3.71318715187768,3.71320203925433,3.71321692663097,3.71323181400762,3.71324670138426,3.71326158876091,3.71327647613755,3.71329136351420,3.71330625089085,3.71332113826749,3.71333602564414,3.71335091302078,3.71336580039743,3.71338068777407,3.71339557515072,3.71341046252736,3.71342534990401,3.71344023728066,3.71345512465730,3.71347001203395,3.71348489941059,3.71349978678724,3.71351467416388,3.71352956154053,3.71354444891717,3.71355933629382,3.71357422367047,3.71358911104711,3.71360399842376,3.71361888580040,3.71363377317705,3.71364866055369,3.71366354793034,3.71367843530699,3.71369332268363,3.71370821006028,3.71372309743692,3.71373798481357,3.71375287219021,3.71376775956686,3.71378264694350,3.71379753432015,3.71381242169680,3.71382730907344,3.71384219645009,3.71385708382673,3.71387197120338,3.71388685858002,3.71390174595667,3.71391663333331,3.71393152070996,3.71394640808661,3.71396129546325,3.71397618283990,3.71399107021654,3.71400595759319,3.71402084496983,3.71403573234648,3.71405061972312,3.71406550709977,3.71408039447642,3.71409528185306,3.71411016922971,3.71412505660635,3.71413994398300,3.71415483135964,3.71416971873629,3.71418460611294,3.71419949348958,3.71421438086623,3.71422926824287,3.71424415561952,3.71425904299616,3.71427393037281,3.71428881774945,3.71430370512610,3.71431859250275,3.71433347987939,3.71434836725604,3.71436325463268,3.71437814200933,3.71439302938597,3.71440791676262,3.71442280413926,3.71443769151591,3.71445257889256,3.71446746626920,3.71448235364585,3.71449724102249,3.71451212839914,3.71452701577578,3.71454190315243,3.71455679052907,3.71457167790572,3.71458656528237,3.71460145265901,3.71461634003566,3.71463122741230,3.71464611478895,3.71466100216559,3.71467588954224,3.71469077691888,3.71470566429553,3.71472055167218,3.71473543904882,3.71475032642547,3.71476521380211,3.71478010117876,3.71479498855540,3.71480987593205,3.71482476330870,3.71483965068534,3.71485453806199,3.71486942543863,3.71488431281528,3.71489920019192,3.71491408756857,3.71492897494521,3.71494386232186,3.71495874969851,3.71497363707515,3.71498852445180,3.71500341182844,3.71501829920509,3.71503318658173,3.71504807395838,3.71506296133502,3.71507784871167,3.71509273608832,3.71510762346496,3.71512251084161,3.71513739821825,3.71515228559490,3.71516717297154,3.71518206034819,3.71519694772483,3.71521183510148,3.71522672247813,3.71524160985477,3.71525649723142,3.71527138460806,3.71528627198471,3.71530115936135,3.71531604673800,3.71533093411465,3.71534582149129,3.71536070886794,3.71537559624458,3.71539048362123,3.71540537099787,3.71542025837452,3.71543514575116,3.71545003312781,3.71546492050446,3.71547980788110,3.71549469525775,3.71550958263439,3.71552447001104,3.71553935738768,3.71555424476433,3.71556913214097,3.71558401951762,3.71559890689427,3.71561379427091,3.71562868164756,3.71564356902420,3.71565845640085,3.71567334377749,3.71568823115414,3.71570311853078,3.71571800590743,3.71573289328408,3.71574778066072,3.71576266803737,3.71577755541401,3.71579244279066,3.71580733016730,3.71582221754395,3.71583710492060,3.71585199229724,3.71586687967389,3.71588176705053,3.71589665442718,3.71591154180382,3.71592642918047,3.71594131655711,3.71595620393376,3.71597109131041,3.71598597868705,3.71600086606370,3.71601575344034,3.71603064081699,3.71604552819363,3.71606041557028,3.71607530294693,3.71609019032357,3.71610507770022,3.71611996507686,3.71613485245351,3.71614973983015,3.71616462720680,3.71617951458344,3.71619440196009,3.71620928933673,3.71622417671338,3.71623906409003,3.71625395146667,3.71626883884332,3.71628372621996,3.71629861359661,3.71631350097325,3.71632838834990,3.71634327572655,3.71635816310319,3.71637305047984,3.71638793785648,3.71640282523313,3.71641771260977,3.71643259998642,3.71644748736306,3.71646237473971,3.71647726211636,3.71649214949300,3.71650703686965,3.71652192424629,3.71653681162294,3.71655169899958,3.71656658637623,3.71658147375288,3.71659636112952,3.71661124850617,3.71662613588281,3.71664102325946,3.71665591063610,3.71667079801275,3.71668568538939,3.71670057276604,3.71671546014268,3.71673034751933,3.71674523489598,3.71676012227262,3.71677500964927,3.71678989702591,3.71680478440256,3.71681967177920,3.71683455915585,3.71684944653250,3.71686433390914,3.71687922128579,3.71689410866243,3.71690899603908,3.71692388341572,3.71693877079237,3.71695365816901,3.71696854554566,3.71698343292231,3.71699832029895,3.71701320767560,3.71702809505224,3.71704298242889,3.71705786980553,3.71707275718218,3.71708764455883,3.71710253193547,3.71711741931212,3.71713230668876,3.71714719406541,3.71716208144205,3.71717696881870,3.71719185619534,3.71720674357199,3.71722163094863,3.71723651832528,3.71725140570193,3.71726629307857,3.71728118045522,3.71729606783186,3.71731095520851,3.71732584258515,3.71734072996180,3.71735561733845,3.71737050471509,3.71738539209174,3.71740027946838,3.71741516684503,3.71743005422167,3.71744494159832,3.71745982897496,3.71747471635161,3.71748960372826,3.71750449110490,3.71751937848155,3.71753426585819,3.71754915323484,3.71756404061148,3.71757892798813,3.71759381536477,3.71760870274142,3.71762359011807,3.71763847749471,3.71765336487136,3.71766825224800,3.71768313962465,3.71769802700129,3.71771291437794,3.71772780175458,3.71774268913123,3.71775757650788,3.71777246388452,3.71778735126117,3.71780223863781,3.71781712601446,3.71783201339110,3.71784690076775,3.71786178814440,3.71787667552104,3.71789156289769,3.71790645027433,3.71792133765098,3.71793622502762,3.71795111240427,3.71796599978091,3.71798088715756,3.71799577453421,3.71801066191085,3.71802554928750,3.71804043666414,3.71805532404079,3.71807021141743,3.71808509879408,3.71809998617072,3.71811487354737];

% 50x50

eps_all = [146.187045931271,146.193018722263,146.198991513255,146.204964304248,146.210937095240,146.216909886232,146.222882677225,146.228855468217,146.234828259209,146.240801050201,146.246773841194,146.252746632186,146.258719423178,146.264692214170,146.270665005163,146.276637796155,146.282610587147,146.288583378139,146.294556169132,146.300528960124,146.306501751116,146.312474542109,146.318447333101,146.324420124093,146.330392915085,146.336365706078,146.342338497070,146.348311288062,146.354284079054,146.360256870047,146.366229661039,146.372202452031,146.378175243024,146.384148034016,146.390120825008,146.396093616000,146.402066406993,146.408039197985,146.414011988977,146.419984779969,146.425957570962,146.431930361954,146.437903152946,146.443875943938,146.449848734931,146.455821525923,146.461794316915,146.467767107908,146.473739898900,146.479712689892];
sig_all = [3.71068607260122,3.71083768078380,3.71098928896637,3.71114089714895,3.71129250533152,3.71144411351409,3.71159572169667,3.71174732987924,3.71189893806182,3.71205054624439,3.71220215442697,3.71235376260954,3.71250537079212,3.71265697897469,3.71280858715726,3.71296019533984,3.71311180352241,3.71326341170499,3.71341501988756,3.71356662807014,3.71371823625271,3.71386984443529,3.71402145261786,3.71417306080043,3.71432466898301,3.71447627716558,3.71462788534816,3.71477949353073,3.71493110171331,3.71508270989588,3.71523431807846,3.71538592626103,3.71553753444360,3.71568914262618,3.71584075080875,3.71599235899133,3.71614396717390,3.71629557535648,3.71644718353905,3.71659879172163,3.71675039990420,3.71690200808677,3.71705361626935,3.71720522445192,3.71735683263450,3.71750844081707,3.71766004899965,3.71781165718222,3.71796326536480,3.71811487354737];

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