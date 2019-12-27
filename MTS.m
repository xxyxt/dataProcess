%% 2��ѹ��Ƭ�������20mm����
clc;
clear;
global Fn f lamda slopAndIntercept sigmaAfterMove hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V_caculate V_plot T_caculate mdata T starPoint sigma_caculate sigma_model endPoint T_model startTime endTime fileName inputFilePath outputFilePath templateSigma
templateSigma='sigmaCaculate' ;       % ģ���ļ�?1?7??
lamda=-4.6343e-10;
beita=3.4088*10^-8;
hp=5*10^-3;                      % ���
rp=10*10^-3;                     % �뾶
d33=650*10^-12;                  % ѹ��ϵ��
varepsilon33=4000*8.854*10^-12;  % ���ϵ��
S0=pi*(rp^2);                    % ѹ��Ƭ�����
P=0.7*10^6;                      % ����Ӧ��
d=10*10^-3;                      % �������
s=pi*(0.15^2);                   % �������

% starPoint=1;                     % ��ȡ������ʼ?1?7??
% endPoint=1451;                   % ��ȡ���ݽ���?1?7??
% inputFilePath='D:\matlabworkbench\MTS\3.xlsx';
% outputFilePath='D:\workBench\matlabworkbench\data\\';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T_caculate,V_caculate,V_plot]=ReadFile();                              % ��ȡʵ������ T_caculateʱ�����  V1��ѹ����
% T=1/f;                           % ��������
% w1=-varepsilon33/(d33*hp);       % ΢�ַ���ϵ��
% w2=-1/(d33*pi*rp^2*R);           % ΢�ַ���ϵ�� 
% % w1=-beita/(hp*lamda);
% % w2=-1/(pi*rp^2*R*lamda);
% T_caculate=T_caculate-T_caculate(1);                     % �����ʼʱ��Ϊ��


% V_caculate=V_caculate-0.2334*0.8;

% % V_caculate=V_caculate-0.2480*2.2;
% sigma_model=PlotCaculateSigma();                         % PlotCaculateSigma ģ��Ӧ���������  ģ��
% sigma_caculate=PlotRealSigma();                          % PlotRealSigma     ����Ӧ���������

% [sigma_simulation,time_simulation]=findMaxValueHalfCycle();
% slopAndIntercept=polyfit(time_simulation,sigma_simulation,1);
% moveDown();
% simulation_time=linspace(min(time_simulation),max(time_simulation));
% simulation_sigma=polyval(slopAndIntercept,simulation_time);

% [startTime,endTime]=find(sigma_model,sigma_caculate);    % ������λ?1?7??
% DealOrignalSigma();
% WriteFile(fileName);
% mdata=[T_caculate(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint),T_caculate(starPoint:endPoint),sigma_model(startTime:endTime)];
% OriginPlot();
%��ͼ
% figure(1)
% plot(T_caculate(starPoint:endPoint),V_caculate,'r.');
% figure(2);                       
% plot(T_caculate(starPoint:endPoint),sigma_model(starPoint:endPoint));
% hold on;
% plot(T_caculate(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
% plot(T_caculate(starPoint:endPoint),sigmaAfterMove(starPoint:endPoint));
% plot(time_simulation,sigma_simulation,'*',simulation_time,simulation_sigma);

% һ�δ�����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='D:\matlabworkbench\MTS\4\';
[ files,num ] = get_all_files( path );
hwait=waitbar(0,'��ȴ�>>>>>>>>');
fList=[1,5,10,8];
FnList=[36000,45000,54000,63000,72000];
j=1;
k=1;
for i=1:1:5
    R0=10*10^6;
    [ T_caculate,V_caculate ]=chooseData(files{i});
    w1=-varepsilon33/(d33*hp);          % ΢�ַ���ϵ��
    w2=-1/(d33*pi*rp^2*R0);             % ΢�ַ���ϵ��  
    fileNameList=split(files{i},'\');
    fileName=[fileNameList{4}];
    T_caculate=T_caculate-T_caculate(1);                     % �����ʼʱ��Ϊ��

    if mod(i,5)==0
        j=j+1;
        k=5;
    end
    if mod(i,5)~=0
        k=mod(i,5);
    end     
    T=1/fList(j);                          % ��������                      
    sigma_model=PlotCaculateSigma(FnList(k),fList(j));       % PlotCaculateSigma ģ��Ӧ���������
    sigma_caculate=PlotRealSigma();                          % PlotRealSigma     ����Ӧ���������
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

    % [startTime,endTime]=find(sigma_model,sigma_caculate);    % ������λ��
    % DealOrignalSigma();                                      % �����ʼֵ����
    % WriteFile(fileName);                                     % ������д���ļ�
    % mdata=[T_caculate,sigma_caculate,T_caculate,sigma_model(startTime:endTime)]; % ��ʼ��origin��ͼ����
    % OriginPlot();                                            % ����origin��ͼ
    
end
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
    % Ѱ�Ҹ�������
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
%���㷴��Ӧ��
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
    % sigma_caculate=-sigma_caculate;
    
    
end
%����ģ��Ӧ��
function sigma_model=PlotCaculateSigma(Fn,f)
    global T_caculate v rp P S0 T_model s
    N=length(T_caculate);
    sigma_model=zeros(N,1);
    for i=1:1:N
        sigma_model(i,1)=-(Fn*sin(2*pi*f*T_caculate(i,1)+pi/2)/s-Fn/s);
    end
end
% %����ģ��Ӧ��
% function sigma_model=PlotCaculateSigma()
%     global T_caculate v rp P S0 T_model Fn f s
%     N=length(T_caculate);
%     sigma_model=zeros(N,1);
%     for i=1:1:N
%         sigma_model(i,1)=-(Fn*sin(2*pi*f*T_caculate(i,1)+pi/2)/s-Fn/s);
%     end
% end
% ƽ��Ӧ��
function moveDown()
    global sigma_caculate sigmaAfterMove T_caculate slopAndIntercept
    len=length(sigma_caculate);
    for i=1:1:len
        sigmaAfterMove(i,1)=sigma_caculate(i,1)-slopAndIntercept(1,1)*T_caculate(i);
    end
end

% �ҳ�T/2�����ڵ�?1?7??��?1?7??1?7��?1?7??
function [sigma_simulation,time_simulation]=findMaxValueHalfCycle()
    global T sigma_caculate T_caculate
    len=length(T_caculate);
    limit=fix(len/T*0.004);
    star=1;
    sigma_simulation=zeros(limit,1);
    time_simulation=zeros(limit,1);
    for i=1:1:limit
        maxNum=0;
        time=0;
        flag=i*T/0.004;
        for j=star:1:flag
            if sigma_caculate(j)>maxNum
                maxNum=sigma_caculate(j);
                time=T_caculate(j);
            end
        end
        star=star+T/0.004;
        sigma_simulation(i,1)=maxNum;
        time_simulation(i,1)=time;
    end
    
end

%������λ?1?7??
function [startTime,endTime]=find(sigma_model,sigma_caculate)
    % 0-138s�ڵ�sigma_ex ��Ӧ��ʱ��t_0 ?1?7??��Ӧ?1?7?? max_sigma_model
    % 0-138s�ڵ�sigma ��Ӧ��ʱ��t_1 ?1?7??��Ӧ?1?7?? max_sigma_caculate
    % ��� diffValue=max_sigma_model-max_sigma_caculate ��?1?7??1?7��С��ʱ�� ��Ӧ��ʱ���?1?7?? distTimeValue=t_0-t_1
    % t>0 ?1?7?? startTime-t*100
    % t<0 ?1?7?? startTime+t*100
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
%��ȡ�ļ�
function [T_caculate,V_caculate,V_plot]=ReadFile()
    global starPoint endPoint R d33 fileName inputFilePath V_plot f Fn
    % ALL0041��1hz-45KN��
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=1;
    % Fn=45*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['F' num2str(starPoint) ':' 'F' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0041��10hz-72KN��
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=10;
    % Fn=72*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['I' num2str(starPoint) ':' 'I' num2str(endPoint)]);
    % V_caculate=V_plot;

    
    % ALL0040��1hz-36KN��
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=1;
    % Fn=36*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0040��1hz-54KN��
    % R=10*10^6
    % endPoint=2249;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=1;
    % Fn=54*10^3;
	% T_caculate=xlsread(inputFilePath,'Sheet1',['K' num2str(starPoint) ':' 'K' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0043��1hz-63KN��
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=1;
    % Fn=63*10^3;
    % T_caculate=xlsread(inputFilePath,'Sheet1',['N' num2str(starPoint) ':' 'N' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['O' num2str(starPoint) ':' 'O' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0044��1hz-72KN��
    % R=10*10^6
    % endPoint=2250;
    % d33=600*10^-12;                  % ѹ��ϵ��
    % varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    % f=1;
    % Fn=72*10^3;
    % T_caculate=xlsread(inputFilePath,'Sheet1',['Q' num2str(starPoint) ':' 'Q' num2str(endPoint)]);
    % V_plot=xlsread(inputFilePath,'Sheet1',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % V_caculate=V_plot;

    % ALL0045��5hz-36KN��
    R=10*10^6
    endPoint=2250;
    d33=600*10^-12;                  % ѹ��ϵ��
    varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    f=5;
    Fn=36*10^3;
    T_caculate=xlsread(inputFilePath,'Sheet1',['T' num2str(starPoint) ':' 'T' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['U' num2str(starPoint) ':' 'U' num2str(endPoint)]);
    V_caculate=V_plot;

    % ALL0045��5hz-45KN��
    R=10*10^6
    endPoint=2250;
    d33=600*10^-12;                  % ѹ��ϵ��
    varepsilon33=4000*8.854*10^-12;  % ���ϵ�� 
    f=5;
    Fn=45*10^3;
    T_caculate=xlsread(inputFilePath,'Sheet1',['W' num2str(starPoint) ':' 'W' num2str(endPoint)]);
    V_plot=xlsread(inputFilePath,'Sheet1',['X' num2str(starPoint) ':' 'X' num2str(endPoint)]);
    V_caculate=V_plot;

end
