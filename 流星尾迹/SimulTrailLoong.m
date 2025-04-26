FileNa = -1;
FileK = -1;
FileFe = -1;
FileNi = -1;
% FileCa = -1;
% FileCap = -1;
FileRay = -1;

NajNum = 1; NaiNum = 1;NaDex = 1/1;
KjNum = 1; KiNum = 1;KDex = 1/2;
FejNum = 1; FeiNum = 1;FeDex = 1/2;
NijNum = 1; NiiNum = 1;NiDex = 1/1;
CajNum = 1; CaiNum = 1;CaDex = 1/1;
CapjNum = 1; CapiNum = 1;CapDex = 1/1;
RayjNum = 1; RayiNum = 1;RayDex = KDex;

% jnK = 559;
% jnFe = 571;
% jnNi = 282;
jnCa= 565;
jnCap = 423;

TimeSum = TimeXCap;
YL = 75;
YR = 115;

PF = jnCap

jt = PF'

TtkUp = 1.5/1440;
TtkDp = 5.5/1440;

Ttk = 1/1440;

HitAlt = 81.56;

start_time = datenum(TimeSum(:,jt))-TtkUp;
end_time = datenum(TimeSum(:,jt))+TtkDp;

% fC = figure('name','Sum','position',[50,100,725,367.5])

% NaV

% subplot(3,3,1); 800 400
figure('name','全队列','position',[10 50 725 367.5]);

if FileNa >= 0
figure('name','Na')
    PhX = PhjSumNaV(KM75NaV:KM115NaV,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaV(KM75NaV:KM115NaV,:).*((1/1440)./PhMax).*NaDex;
    NaV = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaV(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaV,AltitudeNaV(KM75NaV:KM115NaV,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    NaVStr  = "NaV " + RunDate;
    title(NaVStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% NaN

% subplot(3,3,4);

if FileNa >= 0
figure('name','Na2')
    PhX = PhjSumNaM(KM75NaM:KM115NaM,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaM(KM75NaM:KM115NaM,:).*((1/1440)./PhMax).*NaDex;
    NaM = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaM(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaM,AltitudeNaM(KM75NaM:KM115NaM,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    NaMStr  = "NaM " + RunDate;
    title(NaMStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% NaE

% subplot(3,3,7);

if FileNa >= 0
figure('name','Na3')
    PhX = PhjSumNaZ(KM75NaZ:KM115NaZ,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumNaZ(KM75NaZ:KM115NaZ,:).*((1/1440)./PhMax).*NaDex;
    NaZ = PhMin;
    TimeD = datenum(TimeXNa);
    for jn = 1:size(PhMin,2)
        NaZ(:,jn) = datenum(TimeXNa(:,jn)) + PhMin(:,jn);
    end
    plot(NaZ,AltitudeNaZ(KM75NaZ:KM115NaZ,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    NaZStr  = "NaZ " + RunDate;
    title(NaZStr)

    colorX = ['g','b','r'];
    ColorD = repmat(colorX,1,FileNa/3);
    h = findobj(gca, 'Type', 'line');
    for i = 1:numel(h)
        set(h(i), 'Color', ColorD(i));
    end

end

% K

% subplot(3,3,2);

if FileK >= 0
% figure('name','K')
subplot(4,1,1)
    PhX = DenK(KM75K:KM115K,jnK);
    PhMax = max(PhX, [], 1);
    PhMin = DenK(KM75K:KM115K,:).*((1/1440)./PhMax).*KDex;
    K = PhMin;
    TimeD = datenum(TimeXK);
    for jn = 1:size(PhMin,2)
        K(:,jn) = datenum(TimeXK(:,jn)) + PhMin(:,jn);
    end
    
    plot([start_time end_time],[90.56 90.56],'--k','linewidth',0.5)
    hold on
    
    plot(K,AltitudeK(KM75K:KM115K,:),'linewidth',1.5)
    grid on;
    ylabel('Altitude (km)')
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM:ss','keeplimits');
    ax = gca;
    ax.XMinorTick = 'on';
%     ax.XAxis.MinorTickValues = 0:0.1:10;
    KStr  = "K " + RunDate;
    title(KStr)
    set(gca, 'fontsize', 11);
    
end
% MaxK = PhMax;
% disp({'MaxK is ',MaxK})

% subplot(3,3,5);

if FileFe >= 0
% figure('name','Fe')
subplot(4,1,2)
    PhX = DenFe(KM75Fe:KM115Fe,jnFe);
    PhMax = max(PhX, [], 1);
    PhMin = DenFe(KM75Fe:KM115Fe,:).*((1/1440)./PhMax).*FeDex;
    Fe = PhMin;
    TimeD = datenum(TimeXFe);
    for jn = 1:size(PhMin,2)
        Fe(:,jn) = datenum(TimeXFe(:,jn)) + PhMin(:,jn);
    end
    
    plot([start_time end_time],[90.56 90.56],'--k','linewidth',0.5)
    hold on
    
    plot(Fe,AltitudeFe(KM75Fe:KM115Fe,:),'linewidth',1.5)
    grid on;
    ylabel('Altitude (km)')
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM:ss','keeplimits');
    ax = gca;
    ax.XMinorTick = 'on';
    FeStr  = "Fe " + RunDate;
    title(FeStr)
    set(gca, 'fontsize', 11);
    
end
% MaxFe = PhMax;
% disp({'MaxFe is ',MaxFe})

% subplot(3,3,8);

if FileNi >= 0
% figure('name','Ni')
subplot(4,1,3)
    PhX = DenNi(KM75Ni:KM115Ni,jnNi);
    PhMax = max(PhX, [], 1);
    PhMin = DenNi(KM75Ni:KM115Ni,:).*((1/1440)./PhMax).*NiDex;
    Ni = PhMin;
    TimeD = datenum(TimeXNi);
    for jn = 1:size(PhMin,2)
        Ni(:,jn) = datenum(TimeXNi(:,jn)) + PhMin(:,jn);
    end
    
    plot([start_time end_time],[90.56 90.56],'--k','linewidth',0.5)
    hold on
    
    plot(Ni,AltitudeNi(KM75Ni:KM115Ni,:),'linewidth',1.5)
    grid on;
    ylabel('Altitude (km)')
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM:ss','keeplimits');
    ax = gca;
    ax.XMinorTick = 'on';
    NiStr  = "Ni " + RunDate;
    title(NiStr)
    set(gca, 'fontsize', 11);
    
end
% MaxNi = PhMax;
% disp({'MaxNi is ',MaxNi})

% subplot(3,3,3);

if FileCa >= 0
% figure('name','Ca')
subplot(2,1,1)
    PhX = DenCa(KM75Ca:KM115Ca,jnCa);
    PhMax = max(PhX, [], 1);
    PhMin = DenCa(KM75Ca:KM115Ca,:).*((1/1440)./PhMax).*CaDex;
    Ca = PhMin;
    TimeD = datenum(TimeXCa);
    for jn = 1:size(PhMin,2)
        Ca(:,jn) = datenum(TimeXCa(:,jn)) + PhMin(:,jn);
    end
    
    plot([start_time end_time],[HitAlt HitAlt],'--k','linewidth',0.5)
    hold on;
    
    plot(Ca,AltitudeCa(KM75Ca:KM115Ca,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    ylabel('Altitude (km)')
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM:ss','keeplimits');
    ax = gca;
    ax.XMinorTick = 'on';
    CaStr  = "Ca " + RunDate;
    title(CaStr)
    set(gca, 'fontsize', 11);
    
end
MaxCa = PhMax;
disp({'MaxCa is ',MaxCa})

% subplot(3,3,6);

if FileCap >= 0
% figure('name','Cap')
subplot(2,1,2)
    PhX = DenCap(KM75Cap:KM115Cap,jnCap);
    PhMax = max(PhX, [], 1);
    PhMin = DenCap(KM75Cap:KM115Cap,:).*((1/1440)./PhMax).*CapDex;
    Cap = PhMin;
    TimeD = datenum(TimeXCap);
    for jn = 1:size(PhMin,2)
        Cap(:,jn) = datenum(TimeXCap(:,jn)) + PhMin(:,jn);
    end
    
    plot([start_time end_time],[80.24 80.24],'--k','linewidth',0.5)
    hold on;
    
    plot(Cap,AltitudeCap(KM75Cap:KM115Cap,:),'linewidth',1.5)
    grid on;
    ylabel('Altitude (km)')
    xlabel('Time (UT)')
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM:ss','keeplimits');
    ax = gca;
    ax.XMinorTick = 'on';
    CapStr  = "Ca^{+} " + RunDate;
    title(CapStr)
    set(gca, 'fontsize', 11);
    
end
MaxCap = PhMax;
disp({'MaxCap is ',MaxCap})

% subplot(3,3,9);

if FileRay >= 0
figure('name','Ray')
    PhX = PhjSumRay(KM75Ray:KM115Ray,:);
    PhMax = max(PhX, [], 1);
    PhMin = PhjSumRay(KM75Ray:KM115Ray,:).*((1/1440)./PhMax).*RayDex;
    Ray = PhMin;
    TimeD = datenum(TimeXRay);
    for jn = 1:size(PhMin,2)
        Ray(:,jn) = datenum(TimeXRay(:,jn)) + PhMin(:,jn);
    end
    plot(Ray,AltitudeRay(KM75Ray:KM115Ray,:),'linewidth',1.5)
    grid on;
    ylim([YL YR])
    xlim([start_time,end_time]);
    x_ticks = start_time:Ttk:end_time;
    set(gca, 'XTick', x_ticks);
    datetick('x','HH:MM','keeplimits','keepticks');
    RayStr  = "Ray " + RunDate;
    title(RayStr)
    
end

%     % 指定保存路径和文件名
%     figPath = "F:\TrailIMG\SumTrailFind\"+ReadDate+"\";
%     if ~exist(figPath, 'dir')
%         mkdir(figPath);
%     end
%     figV = datestr(TimeSum(:,jt),'yyyymmdd-HH-MM-SS')+"-Fe-"+string(jt);
%     print(fC, '-dpng', '-r300', figPath+figV);
%     close all;



TS = 'Figures ok '