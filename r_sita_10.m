%% 2号压电片埋置深度20mm计算
clc;
clear;
global lamda slopAndIntercept sigmaAfterMove hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V_caculate V_plot T_caculate mdata T starPoint sigma_caculate sigma_model endPoint T_model startTime endTime fileName inputFilePath outputFilePath templateSigma
templateSigma='sigmaCaculate' ;       % 模板文件名
lamda=-4.6343e-10;
beita=3.4088*10^-8;
hp=5*10^-3;                      % 厚度
rp=10*10^-3;                     % 半径
d33=650*10^-12;                  % 压电系数
varepsilon33=3850*8.854*10^-12;  % 介电系数
S0=pi*(rp^2);                    % 压电片的面积
P=0.7*10^6;                      % 荷载应力
d=10*10^-3;                      % 埋置深度
s=50*20*10^-6;                   % 作用面积
a=(s/pi)^(1/2);                  % 等效半径
% alpha=1-(1+(a/d)^2)^(-3/2);      % 衰减系数
alpha=0.90208;
v=200*10^-3;                     % 轮子速度
T_model;                         % 应力模型时间矩阵
T=276;                           % 加载周期
starPoint=1;                     % 读取数据起始点
endPoint=1451;                   % 读取数据结束点
inputFilePath='C:\\Users\\tx\Desktop\\20191111\\data\\OriginalData\\1.xlsx';
outputFilePath='D:\\workBench\\matlabworkbench\\data\\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T_caculate,V_caculate,V_plot]=ReadFile();                              % 读取实验数据 T_caculate时间矩阵  V1电压矩阵
% w1=-varepsilon33/(d33*hp);       % 微分方程系数
% w2=-1/(d33*pi*rp^2*R);           % 微分方程系数 
w1=-beita/(hp*lamda);
w2=-1/(pi*rp^2*R*lamda);
T_caculate=T_caculate-T_caculate(1);                     % 处理初始时间为零


% V_caculate=V_caculate-0.2334*0.8;
% % V_caculate=V_caculate-0.2480*2.2;
sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma 模型应力矩阵计算  模型
sigma_caculate=PlotRealSigma();                          % PlotRealSigma     反演应力矩阵计算

[sigma_simulation,time_simulation]=findMaxValueHalfCycle();
slopAndIntercept=polyfit(time_simulation,sigma_simulation,1);
moveDown();
simulation_time=linspace(min(time_simulation),max(time_simulation));
simulation_sigma=polyval(slopAndIntercept,simulation_time);

[startTime,endTime]=find(sigma_model,sigma_caculate);    % 处理相位差
DealOrignalSigma();
WriteFile(fileName);
mdata=[T_model(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint),T_model(starPoint:endPoint),sigma_model(startTime:endTime)];
OriginPlot();
%绘图
figure(1)
plot(T_model(starPoint:endPoint),V_caculate,'r.');
figure(2);                       
plot(T_model(starPoint:endPoint),sigma_model(startTime:endTime));
hold on;
plot(T_model(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
plot(T_model(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint));
plot(time_simulation,sigma_simulation,'*',simulation_time,simulation_sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OriginPlot()  
    global mdata templateSigma fileName T_caculate V_plot sigmaAfterMove
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
    sigma_caculate=sigma_caculate/10^6;
    
    
end
%计算模型应力
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
        % 远去
        % for i=0.21:0.01:0.74
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.02+v*(i-0.2))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 回来
        % for i=0.75:0.01:1.28
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.128-v*(i-0.74))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 远去
        % for i=1.49:0.01:2.12
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.02+v*(i-1.48))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end
        % % 回来
        % for i=2.13:0.01:2.76
        %     sigma_model(sigma_ex_flag)=(3*0.7*10^6*0.02*0.05*0.01^3)/(2*pi*((0.148-v*(i-2.12))^2+0.01^2)^2.5)/10^6;
        %     sigma_ex_flag=sigma_ex_flag+1;
        %     if sigma_ex_flag>N
        %         break;
        %     end
        % end

    end
end
% 平移应力
function moveDown()
    global sigma_caculate sigmaAfterMove T_caculate slopAndIntercept
    len=length(sigma_caculate);
    for i=1:1:len
        sigmaAfterMove(i,1)=sigma_caculate(i,1)-slopAndIntercept(1,1)*T_caculate(i);
    end
end

% 找出T/2周期内的最大值矩阵
function [sigma_simulation,time_simulation]=findMaxValueHalfCycle()
    global T sigma_caculate T_caculate
    len=length(T_caculate);
    limit=fix(len/(T/2));
    star=1;
    sigma_simulation=zeros(limit,1);
    time_simulation=zeros(limit,1);
    for i=1:1:limit
        maxNum=0;
        time=0;
        flag=i*(T/2);
        for j=star:1:flag
            if sigma_caculate(j)>maxNum
                maxNum=sigma_caculate(j);
                time=T_caculate(j);
            end
        end
        star=star+T/2;
        sigma_simulation(i,1)=maxNum;
        time_simulation(i,1)=time;
    end
    
end

%处理相位差
function [startTime,endTime]=find(sigma_model,sigma_caculate)
    % 0-138s内的sigma_ex 对应的时刻t_0 最大应力 max_sigma_model
    % 0-138s内的sigma 对应的时刻t_1 最大应力 max_sigma_caculate
    % 算出 diffValue=max_sigma_model-max_sigma_caculate 差值最小的时刻 对应的时间差值 distTimeValue=t_0-t_1
    % t>0 则 startTime-t*100
    % t<0 则 startTime+t*100
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
%读取文件
function [T_caculate,V_caculate,V_plot]=ReadFile()
    global starPoint endPoint R d33 fileName inputFilePath V_plot
    % 2-2 滤波12.5 1M  
    % R=1*10^6
    % endPoint=2210;
    % d33=550*10^-12;  
    % fileName='2-2并1M滤波12.5';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;
    
    % 2-2 未滤波 1M
    % R=1*10^6
    % fileName='2-2并1M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 未滤波 2M
    % R=2*10^6
    % d33=450*10^-12; 
    % endPoint=1862;
    % fileName='2-2并2M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 未滤波 4M
    % R=4*10^6
    % d33=750*10^-12; 
    % endPoint=1496;
    % fileName='2-2并4M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['J' num2str(starPoint) ':' 'J' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['K' num2str(starPoint) ':' 'K' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 未滤波 6M
    % R=6*10^6
    % d33=650*10^-12; 
    % endPoint=1203;
    % fileName='2-2并6M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['N' num2str(starPoint) ':' 'N' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 未滤波 8M
    % R=8*10^6
    % d33=550*10^-12; 
    % endPoint=2250;
    % fileName='2-2并8M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['P' num2str(starPoint) ':' 'P' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['Q' num2str(starPoint) ':' 'Q' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 2-2 未滤波 10M
    % R=20*10^6
    % d33=600*10^-12; 
    % endPoint=1396;
    % fileName='2-2串10M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['S' num2str(starPoint) ':' 'S' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['T' num2str(starPoint) ':' 'T' num2str(endPoint)]);
    % V_caculate=2*V_plot;

    % 2-2 未滤波 20M
    % R=30*10^6
    % d33=480*10^-12; 
    % endPoint=1416;
    % fileName='2-2串20M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['V' num2str(starPoint) ':' 'V' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['W' num2str(starPoint) ':' 'W' num2str(endPoint)]);
    % V_caculate=3*V_plot;

    % 2-2 未滤波 30M
    R=40*10^6
    d33=450*10^-12; 
    endPoint=1297;
    fileName='2-2串30M未滤波';
	T_caculate=xlsread(inputFilePath,'Sheet1',['Y' num2str(starPoint) ':' 'Y' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['Z' num2str(starPoint) ':' 'Z' num2str(endPoint)]);
    V_caculate=4*V_plot;

    % 2-2 未滤波 50M
    % R=60*10^6
    % d33=450*10^-12; 
    % endPoint=1408;
    % fileName='2-2串50M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet1',['AB' num2str(starPoint) ':' 'AB' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['AC' num2str(starPoint) ':' 'AC' num2str(endPoint)]);
    % V_caculate=V_plot*6;


	
    % 3-1 滤波12.5 1M
    % R=1*10^6
    % fileName='3-1并1M滤波12.5';
    % endPoint=1968;
    % T_caculate=xlsread(inputFilePath,'Sheet2',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 1M
    % R=1*10^6
    % fileName='3-1并1M未滤波';
    % endPoint=1968;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 2M
    % fileName='3-1并2M未滤波';
    % endPoint=1559;
    % R=2*10^6
	% T_caculate=xlsread(inputFilePath,'Sheet2',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 滤波12.5 2M
    % R=2*10^6
    % endPoint=1559;
    % fileName='3-1并2M滤波12.5';
	% T_caculate=xlsread(inputFilePath,'Sheet2',['I' num2str(starPoint) ':' 'I' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['J' num2str(starPoint) ':' 'J' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 4M
    % R=4*10^6
    % fileName='3-1并4M未滤波';
    % d33=490*10^-12;  
    % endPoint=1914;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 6M
    % R=6*10^6
    % fileName='3-1并6M未滤波';
    % endPoint=1343;
    % d33=450*10^-12;  
	% T_caculate=xlsread(inputFilePath,'Sheet2',['O' num2str(starPoint) ':' 'O' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['P' num2str(starPoint) ':' 'P' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 8M
    % R=8*10^6
    % fileName='3-1并8M未滤波';
    % endPoint=1200;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['S' num2str(starPoint) ':' 'S' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 3-1 未滤波 10M
    % R=20*10^6
    % fileName='3-1串10M未滤波';
    % endPoint=1343;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['U' num2str(starPoint) ':' 'U' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['V' num2str(starPoint) ':' 'V' num2str(endPoint)]);
    % V_caculate=V_plot*2;

    % 3-1 未滤波 20M
    % R=30*10^6
    % d33=500*10^-12; 
    % fileName='3-1串20M未滤波';
    % endPoint=1238;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['X' num2str(starPoint) ':' 'X' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['Y' num2str(starPoint) ':' 'Y' num2str(endPoint)]);
    % V_caculate=V_plot*3;

    % 3-1 未滤波 30M
    % R=40*10^6
    % d33=500*10^-12;
    % fileName='3-1串30M未滤波';
    % endPoint=1463;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['AA' num2str(starPoint) ':' 'AA' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['AB' num2str(starPoint) ':' 'AB' num2str(endPoint)]);
    % V_caculate=V_plot*4;

    % 3-1 未滤波 50M
    % R=60*10^6
    % d33=500*10^-12;
    % fileName='3-1串50M未滤波';
    % endPoint=1250;
	% T_caculate=xlsread(inputFilePath,'Sheet2',['AD' num2str(starPoint) ':' 'AD' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet2',['AE' num2str(starPoint) ':' 'AE' num2str(endPoint)]);
    % V_caculate=V_plot*6;
    
    

    % 4-2 未滤波 1M
    % R=1*10^6
    % fileName='4-2并1M未滤波';
    % endPoint=1451;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 滤波12.5 1M
    % R=1*10^6
    % fileName='4-2并1M滤波12.5';
    % endPoint=1451;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 未滤波 2M
    % R=2*10^6
    % endPoint=1344;
    % fileName='4-2并2M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet3',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 未滤波 4M
    % R=4*10^6
    % endPoint=1442;
    % fileName='4-2并4M未滤波';
	% T_caculate=xlsread(inputFilePath,'Sheet3',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 未滤波 6M
    % R=6*10^6
    % fileName='4-2并6M未滤波';
    % endPoint=1398;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['Q' num2str(starPoint) ':' 'Q' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 未滤波 8M
    % R=8*10^6
    % fileName='4-2并8M未滤波';
    % endPoint=1451;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['V' num2str(starPoint) ':' 'V' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['W' num2str(starPoint) ':' 'W' num2str(endPoint)]);
    % V_caculate=V_plot;

    % 4-2 未滤波 10M
    % R=20*10^6
    % d33=500*10^-12;
    % fileName='4-2串10M未滤波';
    % endPoint=1324;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['Y' num2str(starPoint) ':' 'Y' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['Z' num2str(starPoint) ':' 'Z' num2str(endPoint)]);
    % V_caculate=V_plot*2;

    % 4-2 未滤波 20M
    % R=30*10^6
    % d33=450*10^-12;
    % fileName='4-2串20M未滤波';
    % endPoint=1269;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['AB' num2str(starPoint) ':' 'AB' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['AC' num2str(starPoint) ':' 'AC' num2str(endPoint)]);
    % V_caculate=V_plot*3;

    % 4-2 未滤波 30M
    % R=40*10^6
    % d33=450*10^-12;
    % fileName='4-2串30M未滤波';
    % endPoint=1196;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['AE' num2str(starPoint) ':' 'AE' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['AF' num2str(starPoint) ':' 'AF' num2str(endPoint)]);
    % V_caculate=V_plot*4;

    % 4-2 未滤波 50M
    % R=60*10^6
    % d33=450*10^-12;
    % fileName='4-2串50M未滤波';
    % endPoint=1252;
	% T_caculate=xlsread(inputFilePath,'Sheet3',['AH' num2str(starPoint) ':' 'AH' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet3',['AI' num2str(starPoint) ':' 'AI' num2str(endPoint)]);
    % V_caculate=V_plot*6;
end
