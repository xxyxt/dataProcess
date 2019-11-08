clc;
clear;
global hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V1 T_caculate T starPoint endPoint T_model startTime endTime
hp=5*10-3;                       % �������
rp=10*10^-3;                     % �뾶
d33=650*10^-12;                  % ѹ��ϵ��
varepsilon33=3850*8.854*10^-12;  % ���ϵ��
S0=pi*(rp^2);                    % ѹ��Ƭ�����
P=0.7*10^6;                      % ����Ӧ��
d=10*10^-3;                      % �������
s=50*20*10^-6;                   % �������
a=(s/pi)^(1/2);                  % ��Ч�뾶
alpha=1-(1+(a/d)^2)^(-3/2);      % ˥��ϵ��
v=200*10^-3;                     % �����ٶ�
T_model;                              % Ӧ��ģ��ʱ�����
T=276;                           % ��������
starPoint=1;                     % ��ȡ������ʼ��
endPoint=1000;                   % ��ȡ���ݽ�����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T_caculate,V1]=ReadFile();                              % T_caculateʱ�����  V1��ѹ����                                               
w1=-varepsilon33/(d33*hp);                               % ΢�ַ���ϵ��
w2=-1/(d33*pi*rp^2*R);                                   % ΢�ַ���ϵ��
sigma_model=PlotCaculateSigma();                         % ����Ӧ���������
sigma_caculate=PlotRealSigma();                          % ģ��Ӧ���������
[startTime,endTime]=find(sigma_model,sigma_caculate);    % ������λ��

%��ͼ
figure(1)
plot(T_model(starPoint:endPoint),V1,'r.');
figure(2);                       
plot(T_model(starPoint:endPoint),sigma_model(startTime:endTime));
hold on;
plot(T_model(starPoint:endPoint),sigma_caculate(starPoint:endPoint));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%����ģ��Ӧ��
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
           sigma_model(sigma_ex_flag,1)=(rp^2*acos((rp-x)/rp)-(rp-x)*(x*(2*rp-x))^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.06:0.01:0.1
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*(pi-acos((x-rp)/rp))-(x-rp)*(x*(2*rp-x))^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.11:0.01:0.15
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*(pi-acos((3*rp-x)/rp))-(3*rp-x)*(rp^2-(3*rp-x)^2)^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.16:0.01:0.2
           x=i*v;
           sigma_model(sigma_ex_flag,1)=(rp^2*acos((x-3*rp)/rp)-(x-3*rp)*(rp^2-(x-3*rp)^2)^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
    end
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
            if star==1
                star=star+(T-1);
            else
                star=star+T;
            end
        end
        sigma_caculate(i,1)=temp+w1*V1(i);
    end
    sigma_caculate=-sigma_caculate./alpha;        
    
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
%��ȡ�ļ�
function [T_caculate,V1]=ReadFile()
    global starPoint endPoint R
    % 2-2 �˲�12.5 1M
    R=1*10^6
	T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % 2-2 δ�˲� 1M
    %R=1*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
	
    %3-1 �˲�12.5 1M
    % R=1*10^6
    % T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % 3-1 δ�˲� 1M
    % R=1*10^6
	% T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % 3-1 δ�˲� 2M
    %R=2*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % 3-1 �˲�12.5 2M
    % R=2*10^6
	% T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['I' num2str(starPoint) ':' 'I' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['J' num2str(starPoint) ':' 'J' num2str(endPoint)]);
    % 3-1 δ�˲� 4M
    % R=4*10^6
	% T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % 3-1 δ�˲� 6M
    % R=6*10^6
	% T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['O' num2str(starPoint) ':' 'O' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['P' num2str(starPoint) ':' 'P' num2str(endPoint)]);
    % 3-1 δ�˲� 8M
    % R=8*10^6
	% T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    % V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
	
    % 4-2 δ�˲� 1M
    % R=1*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['B' num2str(starPoint) ':' 'B' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['C' num2str(starPoint) ':' 'C' num2str(endPoint)]);
    % 4-2 �˲�12.5 1M
    % R=1*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['D' num2str(starPoint) ':' 'D' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['E' num2str(starPoint) ':' 'E' num2str(endPoint)]);
    % 4-2 δ�˲� 2M
    % R=2*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['G' num2str(starPoint) ':' 'G' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['H' num2str(starPoint) ':' 'H' num2str(endPoint)]);
    % 4-2 δ�˲� 4M
    % R=4*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['L' num2str(starPoint) ':' 'L' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['M' num2str(starPoint) ':' 'M' num2str(endPoint)]);
    % 4-2 δ�˲� 6M
    % R=6*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['Q' num2str(starPoint) ':' 'Q' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['R' num2str(starPoint) ':' 'R' num2str(endPoint)]);
    % 4-2 δ�˲� 8M
    % R=8*10^6
	%T_caculate=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['V' num2str(starPoint) ':' 'V' num2str(endPoint)]);
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3',['W' num2str(starPoint) ':' 'W' num2str(endPoint)]);
end
