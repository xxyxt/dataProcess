%% 2-4��ѹ��Ƭ�������10mm����
clc;
clear;
global hp rp R d33 varepsilon33 S0 P v s a alpha w1 w2 T    %% ���� 
global sigma_caculate sigma_model  T_model startTime endTime  V1 T_caculate mdata %% ���Ĵ�������
global num outputFilePath template filesNew templatePath storagePath fileName  %% ����ļ�·�����ļ�����
template='sigmaCaculate' ;       % ģ���ļ���
hp=20*10^-3;                     % �������
rp=10*10^-3;                     % �뾶
d33=650*10^-12;                  % ѹ��ϵ��
varepsilon33=3850*8.854*10^-12;  % ���ϵ��
S0=pi*(rp^2);                    % ѹ��Ƭ�����
P=0.7*10^6;                      % ����Ӧ��
s=50*20*10^-6;                   % �������
a=(s/pi)^(1/2);                  % ��Ч�뾶
alpha=1-(1+(a/hp)^2)^(-3/2);      % ˥��ϵ��
v=200*10^-3;                     % �����ٶ�
T=276;                           % ��������
templatePath='D:\workBench\OriginWorkbench\myOrigin';     % originģ��·��
storagePath='D:\workBench\OriginWorkbench\10����2-2';     % originͼƬ�洢·��
path='D:\workBench\���԰�ռ�\10����\2-2\';                % ԭʼ���ݶ�ȡ·��
outputFilePath='D:\workBench\matlabworkbench\10����2-2\'; % origin��ͼ���ݴ洢λ��
R=[1e6,2e6,4e6,6e6,8e6,20e6,30e6,40e6,60e6];              % �������ݵ���ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ files,num ] = get_all_files( path );
filesNew=[files(1:5);files(9:12)];
for i=1:1:9
    [ T_caculate,V1 ]=chooseData(filesNew{i});
    R0=R(i);
    %% �������������
    if i>5
        V1=V1*(R(i)/(1e7));
    end
    w1=-varepsilon33/(d33*hp);          % ΢�ַ���ϵ��
    w2=-1/(d33*pi*rp^2*R0);             % ΢�ַ���ϵ��  
    fileNameList=split(filesNew{i},'\');
    fileName=[fileNameList{5} '_' fileNameList{6}];
    T_caculate=T_caculate-T_caculate(1);                     % �����ʼʱ��Ϊ��
    N=length(T_caculate)+276*3;                              % ��ʼ��ȫ�ֱ���T_model 
    T_model=zeros(N,1);                                     
    sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma ����Ӧ���������
    sigma_caculate=PlotRealSigma();                          % PlotRealSigma     ģ��Ӧ���������
    [startTime,endTime]=find(sigma_model,sigma_caculate);    % ������λ��
    DealOrignalSigma();                                      % �����ʼֵ����
    WriteFile(fileName);                                     % ������д���ļ�
    mdata=[T_caculate,sigma_caculate,T_caculate,sigma_model(startTime:endTime)]; % ��ʼ��origin��ͼ����
    OriginPlot();                                            % ����origin��ͼ
end
% matlab��ͼ
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
%����ģ��Ӧ��
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
%% �����ʼֵ����
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
%���㷴��Ӧ��
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
%������λ��
function [startTime,endTime]=find(sigma_model,sigma_caculate)
    % 0-138s�ڵ�sigma_ex ��Ӧ��ʱ��t_0 ���Ӧ�� max_sigma_model
    % 0-138s�ڵ�sigma ��Ӧ��ʱ��t_1 ���Ӧ�� max_sigma_caculate
    % ��� diffValue=max_sigma_model-max_sigma_caculate ��ֵ��С��ʱ�� ��Ӧ��ʱ���ֵ distTimeValue=t_0-t_1
    % t>0 �� startTime-t*100
    % t<0 �� startTime+t*100
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
% ��ȡpath·����һ��·���ļ��������к�׺Ϊ.csv���ļ�����
% filesΪ���صľ���·������cell sizeΪ���� ��i��Ԫ��Ϊfiles{i}
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
%% ��ԭʼ���ݻ�ú��ʵ���������
%% �������·����path
%% ���� starPoint, endPoint,T_caculate,V1
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