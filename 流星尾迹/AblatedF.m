%% 初始金属成分百分比
MltSi = 13.6;
MltMg = 14.4;
MltFe = 12.1;
MltAl = 1.2;
MltCa = 0.82;
MltNa = 0.8;
MltK = 0.051;
MltTi = 0.0036;

%% 金属成分烧蚀比  37度
AfSi = 0.595;
AfMg = 0.54;
AfFe = 0.713;
AfAl = 0.075;
AfCa = 0.229;
AfNa = 0.924;
AfK = 0.922;
AfTi = 0.228;

%% 金属成分相对原子质量
MsSi = 28.0855;
MsMg = 24.305;
MsFe = 55.845;
MsAl = 26.982;
MsCa = 40.078;
MsNa = 22.99;
MsK = 39.098;
MsTi = 47.867;

%% 基本质量比
BsSi = 10000 * MltSi * AfSi * 0.01;
BsMg = 10000 * MltMg * AfMg * 0.01;
BsFe = 10000 * MltFe * AfFe * 0.01;
BsAl = 10000 * MltAl * AfAl * 0.01;
BsCa = 10000 * MltCa * AfCa * 0.01;
BsNa = 10000 * MltNa * AfNa * 0.01;
BsK = 10000 * MltK * AfK * 0.01;
BsTi = 10000 * MltTi * AfTi * 0.01;

%% 基本数密度比
BmSi = BsSi * MsSi;
BmMg = BsMg * MsMg;
BmFe = BsFe * MsFe;
BmAl = BsAl * MsAl;
BmCa = BsCa * MsCa;
BmNa = BsNa * MsNa;
BmK = BsK * MsK;
BmTi = BsTi * MsTi;
disp(['Si',string(double(BmSi))])
disp(['Mg',string(double(BmMg))])
disp(['Fe',string(double(BmFe))])
disp(['Al',string(double(BmAl))])
disp(['Ca',string(double(BmCa))])
disp(['Na',string(double(BmNa))])
disp(['K',string(double(BmK))])
disp(['Ti',string(double(BmTi))])





































