% 漠河年平均模型温度

% E区
AltE = (90:1:110)';
TemMean = AltE.*0;
for dd = 1:365
    for hh = 1:24
        [TemRay,DenRay] = atmosnrlmsise00(AltE.*1000,53.3,122.7,2023,dd,hh);
        TemMean = TemMean+TemRay(:,2);
    end
end
TE = TemMean./(365*24);
figure('Name','E区年平均温度');
plot(TE,AltE,'-b','linewidth',2);
grid on;
title('Average Annual Temperature (E Layer) in Mohe')
xlabel('Temperature (K)');
ylabel('Altitude (km)');
TESet = mean(TE,1)

% F区
AltF = (200:1:300)';
TemMean = AltF.*0;
for dd = 1:365
    for hh = 1:24
        [TemRay,DenRay] = atmosnrlmsise00(AltF.*1000,53.3,122.7,2023,dd,hh);
        TemMean = TemMean+TemRay(:,2);
    end
end
TF = TemMean./(365*24);
figure('Name','F区年平均温度');
plot(TF,AltF,'-r','linewidth',2);
grid on;
title('Average Annual Temperature (F Layer) in Mohe')
xlabel('Temperature (K)');
ylabel('Altitude (km)');
TFSet = mean(TF,1)




