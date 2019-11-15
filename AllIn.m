%% 2-4号压电片埋置深度10mm计算
clc;
clear;
global hp rp R d33 varepsilon33 S0 P v s a alpha w1 w2 T    %% 参数 
global sigma_caculate sigma_model  T_model startTime endTime  V1 T_caculate mdata %% 核心处理数据
global num outputFilePath template filesNew templatePath storagePath fileName  %% 相关文件路径及文件名称
template='sigmaCaculate' ;       % 模板文件名
hp=20*10^-3;                     % 埋置深度
rp=10*10^-3;                     % 半径
d33=650*10^-12;                  % 压电系数
varepsilon33=3850*8.854*10^-12;  % 介电系数
S0=pi*(rp^2);                    % 压电片的面积
P=0.7*10^6;                      % 荷载应力
s=50*20*10^-6;                   % 作用面积
a=(s/pi)^(1/2);                  % 等效半径
alpha=1-(1+(a/hp)^2)^(-3/2);      % 衰减系数
v=200*10^-3;                     % 轮子速度
T=276;                           % 加载周期
templatePath='D:\workBench\OriginWorkbench\myOrigin';     % origin模板路径
storagePath='D:\workBench\OriginWorkbench\10埋深2-2';     % origin图片存储路径
path='D:\workBench\弹性版空间\10埋深\2-2\';                % 原始数据读取路径
outputFilePath='D:\workBench\matlabworkbench\10埋深2-2\'; % origin绘图数据存储位置
R=[1e6,2e6,4e6,6e6,8e6,20e6,30e6,40e6,60e6];              % 几组数据电阻值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ files,num ] = get_all_files( path );
filesNew=[files(1:5);files(9:12)];
for i=1:1:9
    [ T_caculate,V1 ]=chooseData(filesNew{i});
    R0=R(i);
    %% 处理串联并联情况
    if i>5
        V1=V1*(R(i)/(1e7));
    end
    w1=-varepsilon33/(d33*hp);          % 微分方程系数
    w2=-1/(d33*pi*rp^2*R0);             % 微分方程系数  
    fileNameList=split(filesNew{i},'\');
    fileName=[fileNameList{5} '_' fileNameList{6}];
    T_caculate=T_caculate-T_caculate(1);                     % 处理初始时间为零
    N=length(T_caculate)+276*3;                              % 初始化全局变量T_model 
    T_model=zeros(N,1);                                     
    sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma 反演应力矩阵计算
    sigma_caculate=PlotRealSigma();                          % PlotRealSigma     模型应力矩阵计算
    [startTime,endTime]=find(sigma_model,sigma_caculate);    % 处理相位差
    DealOrignalSigma();                                      % 处理初始值问题
    WriteFile(fileName);                                     % 计算结果写入文件
    mdata=[T_caculate,sigma_caculate,T_caculate,sigma_model(startTime:endTime)]; % 初始化origin绘图数据
    OriginPlot();                                            % 操作origin绘图
end
% matlab绘图
% figure(1)
% plot(T_model(starPoint:endPoint),V1,'r.');
% figure(2);                       
% plot(T_model(starPoint:endPoint),sigma_model(startTime:endTime));
% hold on;
% plot(T_model(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OriginPlot()
    global mdata template fileName storagePath templatePath
    % Obtain Origin COM Server object
    % This will connect to an existing instance of Origin, or create a new one if none exist
    originObj=actxserver('Origin.ApplicationSI');
    invoke(originObj, 'Execute', 'doc -mc 1;');
    % Clear "dirty" flag in Origin to suppress prompt for saving current project
    invoke(originObj, 'IsModified', 'false');
    
    % Load the custom template project
    dir = strcat(templatePath, '\', template, '.opj');
    invoke(originObj, 'Load', dir);
    
    % Send this data over to the Data1 worksheet
    % wks = invoke(originObj, 'FindWorksheet', 'Book1');
    % invoke(wks, 'Name', 'MySheet');
    invoke(originObj, 'PutWorksheet', 'Book1', real(mdata));
    
    % Save graph
    % cmd = 'expGraph type:=emf overwrite := rename tr1.unit := 2 tr1.width := 10000 path:= "';
    cmd = 'expGraph type:=emf overwrite := replace tr1.unit := 2 tr1.width := 10000 path:= "';
    cmd = strcat(cmd, storagePath, '" filename:= "', fileName, '.emf";');
    cmd = strcat(cmd, storagePath, '" filename:= "', fileName, '.emf";');
    invoke(originObj, 'Execute', cmd);
    
    % Release
    release(originObj);
end
%计算模型应力
function sigma_model=PlotCaculateSigma()
	global v rp P S0 T_model T_caculate
    N=length(T_caculate)+276*3;
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
    end
end
%% 处理初始值问题
function DealOrignalSigma()
    global sigma_caculate
    N=length(sigma_caculate(:,1));
    minNum=0;
    for i=1:1:N
        if sigma_caculate(i)<minNum
            minNum=sigma_caculate(i);
        end
    end
    sigma_caculate=sigma_caculate-minNum;
end
%计算反演应力
function sigma_caculate=PlotRealSigma()
	global  w1 w2 T_caculate V1 alpha T
    star=1;
	N=length(T_caculate);
	sigma_caculate=zeros(length(T_caculate),1);
    for i=2:1:N
        temp=0;
        for j=star:1:i-1
            temp=temp+w2*(V1(j+1)+V1(j))*(T_caculate(j+1)-T_caculate(j))/2;
        end
        if mod(i,T)==0
            star=star+(T-1);
        end
        sigma_caculate(i,1)=temp+w1*V1(i);
    end
    sigma_caculate=-sigma_caculate./alpha/10^6;    
    
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
    if VoltageBefore(1)==VoltageBefore(lenBefore)
        targetNum=VoltageBefore(1);
    else
        targetNum=VoltageBefore(lenBefore);
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