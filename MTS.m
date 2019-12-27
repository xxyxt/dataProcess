<<<<<<< HEAD
%% 2号压电片埋置深度20mm计算
clc;
clear;
global Fn f lamda slopAndIntercept sigmaAfterMove hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V_caculate V_plot T_caculate mdata T starPoint sigma_caculate sigma_model endPoint T_model startTime endTime fileName inputFilePath outputFilePath templateSigma
templateSigma='sigmaCaculate' ;       % 模板文件?1?7??
lamda=-4.6343e-10;
beita=3.4088*10^-8;
hp=5*10^-3;                      % 厚度
rp=10*10^-3;                     % 半径
d33=650*10^-12;                  % 压电系数
varepsilon33=4000*8.854*10^-12;  % 介电系数
S0=pi*(rp^2);                    % 压电片的面积
P=0.7*10^6;                      % 荷载应力
d=10*10^-3;                      % 埋置深度
s=pi*(0.15^2);                   % 作用面积

% starPoint=1;                     % 读取数据起始?1?7??
% endPoint=1451;                   % 读取数据结束?1?7??
% inputFilePath='D:\matlabworkbench\MTS\3.xlsx';
% outputFilePath='D:\workBench\matlabworkbench\data\\';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T_caculate,V_caculate,V_plot]=ReadFile();                              % 读取实验数据 T_caculate时间矩阵  V1电压矩阵
% T=1/f;                           % 加载周期
% w1=-varepsilon33/(d33*hp);       % 微分方程系数
% w2=-1/(d33*pi*rp^2*R);           % 微分方程系数 
% % w1=-beita/(hp*lamda);
% % w2=-1/(pi*rp^2*R*lamda);
% T_caculate=T_caculate-T_caculate(1);                     % 处理初始时间为零


% V_caculate=V_caculate-0.2334*0.8;

% % V_caculate=V_caculate-0.2480*2.2;
% sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma 模型应力矩阵计算  模型
% sigma_caculate=PlotRealSigma();                          % PlotRealSigma     反演应力矩阵计算

% [sigma_simulation,time_simulation]=findMaxValueHalfCycle();
% slopAndIntercept=polyfit(time_simulation,sigma_simulation,1);
% moveDown();
% simulation_time=linspace(min(time_simulation),max(time_simulation));
% simulation_sigma=polyval(slopAndIntercept,simulation_time);

% [startTime,endTime]=find(sigma_model,sigma_caculate);    % 处理相位?1?7??
% DealOrignalSigma();
% WriteFile(fileName);
% mdata=[T_caculate(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint),T_caculate(starPoint:endPoint),sigma_model(startTime:endTime)];
% OriginPlot();
%绘图
% figure(1)
% plot(T_caculate(starPoint:endPoint),V_caculate,'r.');
% figure(2);                       
% plot(T_caculate(starPoint:endPoint),sigma_model(starPoint:endPoint));
% hold on;
% plot(T_caculate(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
% plot(T_caculate(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint));
% plot(time_simulation,sigma_simulation,'*',simulation_time,simulation_sigma);

% 一次处理所有数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='D:\matlabworkbench\MTS\4\';
[ files,num ] = get_all_files( path );
hwait=waitbar(0,'请等待>>>>>>>>');
fList=[1,5,10,8];
FnList=[36000,45000,54000,63000,72000];
j=1;
k=1;
for i=1:1:5
    R0=10*10^6;
    [ T_caculate,V_caculate ]=chooseData(files{i});
    w1=-varepsilon33/(d33*hp);          % 微分方程系数
    w2=-1/(d33*pi*rp^2*R0);             % 微分方程系数  
    fileNameList=split(files{i},'\');
    fileName=[fileNameList{4}];
    T_caculate=T_caculate-T_caculate(1);                     % 处理初始时间为零

    if mod(i,5)==0
        j=j+1;
        k=5;
    end
    if mod(i,5)~=0
        k=mod(i,5);
    end     
    T=1/fList(j);                          % 加载周期                      
    sigma_model=PlotCaculateSigma(FnList(k),fList(j));       % PlotCaculateSigma 模型应力矩阵计算
    sigma_caculate=PlotRealSigma();                          % PlotRealSigma     反演应力矩阵计算
    [sigma_simulation,time_simulation]=findMaxValueHalfCycle();
    slopAndIntercept=polyfit(time_simulation,sigma_simulation,1);
    moveDown();
    simulation_time=linspace(min(time_simulation),max(time_simulation));
    simulation_sigma=polyval(slopAndIntercept,simulation_time);
    figure(i);
    plot(T_caculate,V_caculate,'r.');
    figure(21+i);
    hold on;
    plot(T_caculate,sigma_model);
    plot(T_caculate,sigmaAfterMove);

    % [startTime,endTime]=find(sigma_model,sigma_caculate);    % 处理相位差
    % DealOrignalSigma();                                      % 处理初始值问题
    % WriteFile(fileName);                                     % 计算结果写入文件
    % mdata=[T_caculate,sigma_caculate,T_caculate,sigma_model(startTime:endTime)]; % 初始化origin绘图数据
    % OriginPlot();                                            % 操作origin绘图
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OriginPlot()  
    global mdata templateSigma fileName T_caculate V_plot sigmaAfterMove
=======
%% 2鍙峰帇鐢电墖鍩嬬疆娣卞害20mm璁＄畻
clc;
clear;
global slopAndIntercept sigmaAfterMove hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V_caculate V_plot T_caculate mdata T starPoint sigma_caculate sigma_model endPoint T_model startTime endTime fileName inputFilePath outputFilePath templateSigma
templateSigma='sigmaCaculate' ;       % 妯℃澘鏂囦欢鍚�
hp=5*10^-3;                      % 鍘氬害
rp=10*10^-3;                     % 鍗婂緞
d33=650*10^-12;                  % 鍘嬬數绯绘暟
varepsilon33=3850*8.854*10^-12;  % 浠嬬數绯绘暟
S0=pi*(rp^2);                    % 鍘嬬數鐗囩殑闈㈢Н
P=0.7*10^6;                      % 鑽疯浇搴斿姏
d=20*10^-3;                      % 鍩嬬疆娣卞害
% s=50*20*10^-6;                   % 浣滅敤闈㈢Н
% a=(s/pi)^(1/2);                  % 绛夋晥鍗婂緞
s=pi*(150*10^-3)^2;
% alpha=1-(1+(a/d)^2)^(-3/2);      % 琛板噺绯绘暟
% alpha=0.588695387;               % 绾帇妯″瀷
alpha=1;
v=200*10^-3;                     % 杞瓙閫熷害
T_model;                         % 搴斿姏妯″瀷鏃堕棿鐭╅樀
T=276;                           % 鍔犺浇鍛ㄦ湡
starPoint=1;                     % 璇诲彇鏁版嵁璧峰鐐�
endPoint=1451;                   % 璇诲彇鏁版嵁缁撴潫鐐�
inputFilePath='C:\\Users\\tx\Desktop\\20191111\\data\\OriginalData\\3.xlsx';
outputFilePath='D:\\workBench\\matlabworkbench\\data3_20\\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T_caculate,V_caculate,V_plot]=ReadFile();                              % 璇诲彇瀹為獙鏁版嵁 T_caculate鏃堕棿鐭╅樀  V1鐢靛帇鐭╅樀
w1=-varepsilon33/(d33*hp);       % 寰垎鏂圭▼绯绘暟
w2=-1/(d33*pi*rp^2*R);           % 寰垎鏂圭▼绯绘暟  
T_caculate=T_caculate-T_caculate(1);                     % 澶勭悊鍒濆鏃堕棿涓洪浂
% V_caculate=V_caculate-0.0350*0.9;

% sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma 妯″瀷搴斿姏鐭╅樀璁＄畻  妯″瀷
sigma_caculate=PlotRealSigma();                          % PlotRealSigma     鍙嶆紨搴斿姏鐭╅樀璁＄畻
% [sigma_simulation,time_simulation]=findMaxValueHalfCycle();

% slopAndIntercept=polyfit(time_simulation,sigma_simulation,1);
% simulation_time=linspace(min(time_simulation),max(time_simulation));
% simulation_sigma=polyval(slopAndIntercept,simulation_time);
% sigmaAfterMove=moveDown();
% [startTime,endTime]=find(sigma_model,sigma_caculate);    % 澶勭悊鐩镐綅宸�
% DealOrignalSigma();
% WriteFile(fileName);
% mdata=[T_model(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint),T_model(starPoint:endPoint),sigma_model(startTime:endTime)];
% OriginPlot();
%缁樺浘
figure(1)
plot(T_caculate(starPoint:endPoint),V_caculate,'r.');
figure(2);                       
% plot(T_model(starPoint:endPoint),sigma_model(startTime:endTime));
hold on;
plot(T_caculate(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
% plot(T_model(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint));
% plot(time_simulation,sigma_simulation,'*',simulation_time,simulation_sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OriginPlot()  
    global mdata templateSigma fileName T_caculate V_plot
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
    % Obtain Origin COM Server object
    % This will connect to an existing instance of Origin, or create a new one if none exist
    originObj=actxserver('Origin.ApplicationSI');
    invoke(originObj, 'Execute', 'doc -mc 1;');
    % Clear "dirty" flag in Origin to suppress prompt for saving current project
    invoke(originObj, 'IsModified', 'false');
    
    % Load the custom templateSigma project
    fdir= 'D:\workBench\OriginWorkbench\figure20';
    dirSigma = 'D:\workBench\OriginWorkbench\myOrigin';
    dirSigma = strcat(dirSigma, '\', templateSigma, '.opj');
    invoke(originObj, 'Load', dirSigma);
    
    % Send this data over to the Data1 worksheet
    % wks = invoke(originObj, 'FindWorksheet', 'Book1');
    % invoke(wks, 'Name', 'MySheet');
    invoke(originObj, 'PutWorksheet', '[Book1]Sheet1', real(mdata));

    invoke(originObj, 'PutWorksheet', '[Book1]Sheet2', [T_caculate,V_plot]);
    % Save graph
    % cmd = 'expGraph type:=emf overwrite := rename tr1.unit := 2 tr1.width := 10000 path:= "';
    cmd = 'expGraph type:=emf overwrite := replace tr1.unit := 2 tr1.width := 10000 path:= "';
    cmd = strcat(cmd, fdir, '" filename:= "', fileName, '.emf";');
    cmd = strcat(cmd, fdir, '" filename:= "', fileName, '.emf";');
    invoke(originObj, 'Execute', cmd);
    
    % Release
    release(originObj);
end
<<<<<<< HEAD
% 获取path路径下一级路径文件夹下所有后缀为.csv的文件名称
% files为返回的绝对路径名称cell size为个数 第i个元素为files{i}
function [ files,num ] = get_all_files( path )
    listing=dir(path);
    files={};
    for i=1:length(listing)
        f=listing(i);
        if ~strcmp(f.name,'.') && ~strcmp(f.name,'..')
            fileFolder=fullfile(f.name);
            dirOutput=dir(fullfile([path fileFolder],'*.CSV'));
            files=[files;[path fileFolder '\' dirOutput.name]];
        end
    end
    size0=size(files);
    num=size0(1);
end
%% 从原始数据获得合适的数据区间
%% 输入绝对路径：path
%% 返回 starPoint, endPoint,T_caculate,V1
function [ T_caculate,V1 ]=chooseData(fileName)
    data=xlsread(fileName);
    TimeBefore=data(:,3);
    VoltageBefore=data(:,4);
    lenBefore=length(VoltageBefore);
    starPoint=1;
    endPoint=lenBefore;
    % 寻找干扰数据
    maxNum=-1e6;
    minNum=1e6;
    targetNum=-1e6;
    for i=1:1:lenBefore
        if VoltageBefore(i)>maxNum
            maxNum=VoltageBefore(i);
        end
        if VoltageBefore(i)<minNum
            minNum=VoltageBefore(i);
        end
    end
    if VoltageBefore(1)==minNum || VoltageBefore(lenBefore)==minNum
        targetNum=minNum;
    end
    if VoltageBefore(1)==maxNum || VoltageBefore(lenBefore)==maxNum
        targetNum=maxNum;
    end
    for i=1:1:lenBefore
        if VoltageBefore(i) ~= targetNum
            starPoint=i;
            break;
        end
    end
    for i=lenBefore:-1:starPoint
        if VoltageBefore(i) ~= targetNum
            endPoint=i;
            break;
        end
    end
    T_caculate=TimeBefore(starPoint:endPoint);
    V1=VoltageBefore(starPoint:endPoint);
end
function DealOrignalSigma()
    global  sigmaAfterMove
    N=length(sigmaAfterMove(:,1));
    minNum=0;
    for i=1:1:N
        if sigmaAfterMove(i)<minNum
            minNum=sigmaAfterMove(i);
        end
    end
    sigmaAfterMove=sigmaAfterMove-minNum;
end
%计算反演应力
function sigma_caculate=PlotRealSigma()
	global  w1 w2 T_caculate V_caculate
=======

%璁＄畻鍙嶆紨搴斿姏
function sigma_caculate=PlotRealSigma()
	global  w1 w2 T_caculate V_caculate alpha
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
    star=1;
	N=length(T_caculate);
	sigma_caculate=zeros(length(T_caculate),1);
    for i=2:1:N
        temp=0;
        for j=star:1:i-1
            temp=temp+w2*(V_caculate(j+1)+V_caculate(j))*(T_caculate(j+1)-T_caculate(j))/2;
        end
        % if mod(i,T/2)==0
        %     star=star+(T/2-1);
        % end
        sigma_caculate(i,1)=temp+w1*V_caculate(i);
    end
<<<<<<< HEAD
    % sigma_caculate=-sigma_caculate;
    
    
end
%计算模型应力
function sigma_model=PlotCaculateSigma(Fn,f)
    global T_caculate v rp P S0 T_model s
    N=length(T_caculate);
    sigma_model=zeros(N,1);
    for i=1:1:N
        sigma_model(i,1)=-(Fn*sin(2*pi*f*T_caculate(i,1)+pi/2)/s-Fn/s);
    end
end
% %计算模型应力
% function sigma_model=PlotCaculateSigma()
%     global T_caculate v rp P S0 T_model Fn f s
%     N=length(T_caculate);
%     sigma_model=zeros(N,1);
%     for i=1:1:N
%         sigma_model(i,1)=-(Fn*sin(2*pi*f*T_caculate(i,1)+pi/2)/s-Fn/s);
%     end
% end
% 平移应力
function moveDown()
    global sigma_caculate sigmaAfterMove T_caculate slopAndIntercept
=======
    sigma_caculate=-sigma_caculate;
    
end
%璁＄畻妯″瀷搴斿姏
function sigma_model=PlotCaculateSigma()
	global T_caculate v rp P S0 T_model
    N=length(T_caculate)+276*3;
    T_model=zeros(N,1);
	limit=fix((N-64)/138);
    sigma_model=zeros(N,1);
    index=1;
    for i=0:0.01:(N-1)/100
        T_model(index,1)=i;
        index=index+1;
    end
    for j=0:1:limit
        sigma_ex_flag=j*138+64;   
        for i=0.01:0.01:0.05
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*acos((rp-x)/rp)-(rp-x)*(x*(2*rp-x))^(1/2))*P/S0/10^6;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.06:0.01:0.1
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*(pi-acos((x-rp)/rp))-(x-rp)*(x*(2*rp-x))^(1/2))*P/S0/10^6;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.11:0.01:0.15
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*(pi-acos((3*rp-x)/rp))-(3*rp-x)*(rp^2-(3*rp-x)^2)^(1/2))*P/S0/10^6;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.16:0.01:0.2
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*acos((x-3*rp)/rp)-(x-3*rp)*(rp^2-(x-3*rp)^2)^(1/2))*P/S0/10^6;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        %  % 杩滃幓
        %  for i=0.21:0.01:0.74
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.02+v*(i-0.2))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 鍥炴潵
        % for i=0.75:0.01:1.28
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.128-v*(i-0.74))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 杩滃幓
        % for i=1.49:0.01:2.12
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.02+v*(i-1.48))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 鍥炴潵
        % for i=2.13:0.01:2.76
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.148-v*(i-2.12))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
    end
end
function DealOrignalSigma()
    global sigmaAfterMove
    N=length(sigmaAfterMove(:,1));
    minNum=0;
    for i=1:1:N
        if sigmaAfterMove(i)<minNum
            minNum=sigmaAfterMove(i);
        end
    end
    sigmaAfterMove=sigmaAfterMove-minNum;
end
% 骞崇Щ搴斿姏
function sigmaAfterMove=moveDown()
    global sigma_caculate T_caculate slopAndIntercept
    sigmaAfterMove=zeros(length(sigma_caculate),1);
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
    len=length(sigma_caculate);
    for i=1:1:len
        sigmaAfterMove(i,1)=sigma_caculate(i,1)-slopAndIntercept(1,1)*T_caculate(i);
    end
end

<<<<<<< HEAD
% 找出T/2周期内的?1?7??大?1?7??1?7矩?1?7??
function [sigma_simulation,time_simulation]=findMaxValueHalfCycle()
    global T sigma_caculate T_caculate
    len=length(T_caculate);
    limit=fix(len/T*0.004);
=======
% 鎵惧嚭T/2鍛ㄦ湡鍐呯殑鏈€澶у€肩煩闃�
function [sigma_simulation,time_simulation]=findMaxValueHalfCycle()
    global T sigma_caculate T_caculate
    len=length(T_caculate);
    limit=fix(len/(T/2));
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
    star=1;
    sigma_simulation=zeros(limit,1);
    time_simulation=zeros(limit,1);
    for i=1:1:limit
        maxNum=0;
        time=0;
<<<<<<< HEAD
        flag=i*T/0.004;
        for j=star:1:flag
            if sigma_caculate(j)>maxNum
                maxNum=sigma_caculate(j);
                time=T_caculate(j);
            end
        end
        star=star+T/0.004;
=======
        flag=i*(T/2);
        for j=star:1:flag
            if sigma_caculate(j,1)>maxNum
                maxNum=sigma_caculate(j,1);
                time=T_caculate(j,1);
            end
        end
        star=star+(T/2);
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
        sigma_simulation(i,1)=maxNum;
        time_simulation(i,1)=time;
    end
    
end

<<<<<<< HEAD
%处理相位?1?7??
function [startTime,endTime]=find(sigma_model,sigma_caculate)
    % 0-138s内的sigma_ex 对应的时刻t_0 ?1?7??大应?1?7?? max_sigma_model
    % 0-138s内的sigma 对应的时刻t_1 ?1?7??大应?1?7?? max_sigma_caculate
    % 算出 diffValue=max_sigma_model-max_sigma_caculate 差?1?7??1?7最小的时刻 对应的时间差?1?7?? distTimeValue=t_0-t_1
    % t>0 ?1?7?? startTime-t*100
    % t<0 ?1?7?? startTime+t*100
=======
%澶勭悊鐩镐綅宸�
function [startTime,endTime]=find(sigma_model,sigma_caculate)
    % 0-138s鍐呯殑sigma_ex 瀵瑰簲鐨勬椂鍒籺_0 鏈€澶у簲鍔� max_sigma_model
    % 0-138s鍐呯殑sigma 瀵瑰簲鐨勬椂鍒籺_1 鏈€澶у簲鍔� max_sigma_caculate
    % 绠楀嚭 diffValue=max_sigma_model-max_sigma_caculate 宸€兼渶灏忕殑鏃跺埢 瀵瑰簲鐨勬椂闂村樊鍊� distTimeValue=t_0-t_1
    % t>0 鍒� startTime-t*100
    % t<0 鍒� startTime+t*100
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
    global  T_model T T_caculate
    max_sigma_model=0;
    max_sigma_caculate=0;
    t_0=0;
    t_1=0;
    for t=T+1:1:1.5*T
        if sigma_model(t)>max_sigma_model
            max_sigma_model=sigma_model(t);
            t_0=T_model(t);
        end
    end
    for t=1:1:T/2
        if sigma_caculate(t)>max_sigma_caculate
            max_sigma_caculate=sigma_caculate(t);
            t_1=T_caculate(t);
        end
    end
    distTimeValue=t_0-t_1;
    if distTimeValue>0
        startTime=T+1+distTimeValue*100;
        endTime=startTime+length(sigma_caculate)-1;
    elseif distTimeValue<=0
        startTime=T+1-distTimeValue*100;
        endTime=startTime+length(sigma_caculate)-1;
    end
end
function WriteFile(fileName)
    global T_caculate T_model sigma_model sigma_caculate startTime endTime starPoint endPoint outputFilePath
    delete([outputFilePath fileName '.xlsx']);
    xlswrite([outputFilePath fileName '.xlsx'],T_caculate,['A' num2str(starPoint) ':' 'A' num2str(endPoint)] );
    xlswrite([outputFilePath fileName '.xlsx'],sigma_caculate,['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    xlswrite([outputFilePath fileName '.xlsx'],T_model,['C' num2str(starPoint) ':' 'C' num2str(endPoint)] );
    xlswrite([outputFilePath fileName '.xlsx'],sigma_model(startTime:endTime),['D' num2str(starPoint) ':' 'D' num2str(endPoint)] );
end
<<<<<<< HEAD
%读取文件
function [T_caculate,V_caculate,V_plot]=ReadFile()
    global starPoint endPoint R d33 fileName inputFilePath V_plot f Fn
    % ALL0041【1hz-45KN】
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=1;
    % Fn=45*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['F' num2str(starPoint) ':' 'F' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0041【10hz-72KN】
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=10;
    % Fn=72*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['I' num2str(starPoint) ':' 'I' num2str(endPoint)]);
    % V_caculate=V_plot;

    
    % ALL0040【1hz-36KN】
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=1;
    % Fn=36*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0040【1hz-54KN】
    % R=10*10^6
    % endPoint=2249;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=1;
    % Fn=54*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['K' num2str(starPoint) ':' 'K' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0043【1hz-63KN】
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=1;
    % Fn=63*10^3;
    % T_caculate=xlsread(inputFilePath,'Sheet1',['N' num2str(starPoint) ':' 'N' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['O' num2str(starPoint) ':' 'O' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0044【1hz-72KN】
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % 压电系数
    % varepsilon33=4000*8.854*10^-12;  % 介电系数 
    % f=1;
    % Fn=72*10^3;
    % T_caculate=xlsread(inputFilePath,'Sheet1',['Q' num2str(starPoint) ':' 'Q' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0045【5hz-36KN】
    R=10*10^6
    endPoint=2250;
    d33=600*10^-12;                  % 压电系数
    varepsilon33=4000*8.854*10^-12;  % 介电系数 
    f=5;
    Fn=36*10^3;
    T_caculate=xlsread(inputFilePath,'Sheet1',['T' num2str(starPoint) ':' 'T' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['U' num2str(starPoint) ':' 'U' num2str(endPoint)]);
    V_caculate=V_plot;

    % ALL0045【5hz-45KN】
    R=10*10^6
    endPoint=2250;
    d33=600*10^-12;                  % 压电系数
    varepsilon33=4000*8.854*10^-12;  % 介电系数 
    f=5;
    Fn=45*10^3;
    T_caculate=xlsread(inputFilePath,'Sheet1',['W' num2str(starPoint) ':' 'W' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['X' num2str(starPoint) ':' 'X' num2str(endPoint)]);
    V_caculate=V_plot;

=======
%璇诲彇鏂囦欢
function [T_caculate,V_caculate,V_plot]=ReadFile()
    global starPoint endPoint R d33 fileName inputFilePath V_plot
    R=10*10^6
    endPoint=2250;
    d33=650*10^-12;  
    fileName='1';
	T_caculate=xlsread(inputFilePath,'Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['F' num2str(starPoint) ':' 'F' num2str(endPoint)]);
    V_caculate=V_plot;

    % 2-2 鏈护娉� 1M  
    % R=1*10^6
    % endPoint=1456;
    % d33=550*10^-12;  
    % fileName='2-2骞�1M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 婊ゆ尝15 1M  
    % R=1*10^6
    % endPoint=1456;
    % d33=650*10^-12;  
    % fileName='2-2骞�1M婊ゆ尝15';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 鏈护娉� 2M  
    % R=2*10^6
    % endPoint=1875;
    % d33=650*10^-12;  
    % fileName='2-2骞�2M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 婊ゆ尝15 2M  
    % R=2*10^6
    % endPoint=1456;
    % d33=650*10^-12;  
    % fileName='2-2骞�2M婊ゆ尝15';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['I' num2str(starPoint) ':' 'I' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['J' num2str(starPoint) ':' 'J' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 鏈护娉� 4M  
    % R=4*10^6
    % endPoint=2068;
    % d33=450*10^-12;  
    % fileName='2-2骞�4M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['S' num2str(starPoint) ':' 'S' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 鏈护娉� 6M  
    % R=6*10^6
    % endPoint=2048;
    % d33=650*10^-12;  
    % fileName='2-2骞�6M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['O' num2str(starPoint) ':' 'O' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['P' num2str(starPoint) ':' 'P' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 鏈护娉� 8M  
    % R=8*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2骞�8M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 鏈护娉� 3M  
    % R=13*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2涓�3M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AA' num2str(starPoint) ':' 'AA' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AB' num2str(starPoint) ':' 'AB' num2str(endPoint)]);
    % V_caculate=V_plot*1.3;

    % 2-2 鏈护娉� 6M  
    % R=16*10^6
    % endPoint=2249;
    % d33=700*10^-12;  
    % fileName='2-2涓�3M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AD' num2str(starPoint) ':' 'AD' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AE' num2str(starPoint) ':' 'AE' num2str(endPoint)]);
    % V_caculate=V_plot*1.6;

    
    % 2-2 鏈护娉� 10M  
    % R=20*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2涓�10M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AG' num2str(starPoint) ':' 'AG' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AH' num2str(starPoint) ':' 'AH' num2str(endPoint)]);
    % V_caculate=V_plot*2;

    % 2-2 鏈护娉� 20M  
    % R=30*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2涓�20M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AJ' num2str(starPoint) ':' 'AJ' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AK' num2str(starPoint) ':' 'AK' num2str(endPoint)]);
    % V_caculate=V_plot*3;

    % 2-2 鏈护娉� 30M  
    % R=40*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2涓�30M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AM' num2str(starPoint) ':' 'AM' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AN' num2str(starPoint) ':' 'AN' num2str(endPoint)]);
    % V_caculate=V_plot*4;

    % 2-2 鏈护娉� 50M  
    % R=60*10^6
    % endPoint=1356;
    % d33=650*10^-12;  
    % fileName='2-2涓�50M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AP' num2str(starPoint) ':' 'AP' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AQ' num2str(starPoint) ':' 'AQ' num2str(endPoint)]);
    % V_caculate=V_plot*6;

    % 2-2 鏈护娉� 60M  
    % R=70*10^6
    % endPoint=2250;
    % d33=650*10^-12;  
    % fileName='2-2涓�60M鏈护娉�';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AM' num2str(starPoint) ':' 'AM' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AN' num2str(starPoint) ':' 'AN' num2str(endPoint)]);
    % V_caculate=V_plot*7;
>>>>>>> 240a47067f701a122fa660ff29324162f8ea9be5
end
