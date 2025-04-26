XlsName ='F:\TrailIMG\CaCapRatioXie.xlsx';
Tablex = readtable(XlsName,'sheet','Calcu');
TabDatx = table2array(Tablex);
TimeRo = TabDatx(:,1);
CaDenx = TabDatx(:,2);
CapDenx = TabDatx(:,3);
Ratiox = CaDenx./CapDenx;

% cftool(TimeRo,Ratiox);
Fpoly1 = fit(TimeRo,Ratiox, 'poly1');
Coe = coeffvalues(Fpoly1);
Xie = Coe(1)
scatter(TimeRo,Ratiox,50,'b','linewidth',2);
hold on; box on; grid on;
plot(Fpoly1,'-r')
xlabel('Time in seconds');
ylabel('Ratio of Ca vs Ca ion')
Titlestr = 'Case12 | ' + string(Xie);
title(Titlestr);
legend('Ratio','Fit curve')

%% Ratio
XlsName ='F:\TrailIMG\CaCapRatioXie.xlsx';
Tabley = readtable(XlsName,'sheet','Xie');
TabDaty = table2array(Tabley);
Alty = TabDaty(:,1);
Ratioy = TabDaty(:,2);

% cftool(Alty,Ratioy)
flx = fit(Alty,Ratioy,'poly1');
flmix = -0.002058.*Alty+0.276;

scatter(Ratioy,Alty,50,'b','linewidth',2);
% cftool(Alty,Ratioy);
box on; grid on; hold on;
plot(flmix,Alty,'-r','linewidth',2);
hold on;
plot([0 0],[80 100],'-k','linewidth',2);
xlabel('Ratio');
ylabel('Altitude (km)');
title('Slope of Ca/Ca ion ratio variation')
legend('Slope','Fit Curve')









