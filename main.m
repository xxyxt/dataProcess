clc;
clear;
global hp rp R d33 varepsilon33 S0 P v d s a alpha w1 w2 V1 T1 T
hp=5*10-3;                       %埋置深度
rp=10*10^-3;                     %半径
R=8*10^6;                        %电阻
d33=650*10^-12;                  %压电系数
varepsilon33=3850*8.854*10^-12;  %介电系数
S0=pi*(rp^2);                    %压电片的面积
P=0.7*10^6;                      %荷载应力
d=10*10^-3;                      %埋置深度
s=50*20*10^-6;                   %作用面积
a=(s/pi)^(1/2);                  %等效半径
alpha=1-(1+(a/d)^2)^(-3/2);      %衰减系数
v=200*10^-3;                     %轮子速度
w1=-varepsilon33/(d33*hp);       %微分方程系数
w2=-1/(d33*pi*rp^2*R);           %微分方程系数
T=276;                           %周期

[T1,V1]=ReadFile();
sigma_ex=PlotCaculateSigma();    %反演应力矩阵计算
sigma=PlotRealSigma();           %模型应力矩阵计算
figure(1)
plot(T1,V1,'r.');
figure(2);                       %绘图
plot(T1,sigma_ex);
hold on;
plot(T1,sigma);

%计算反演应力
function sigma_ex=PlotCaculateSigma()
	global T1 v rp P S0
	N=length(T1);
	limit=fix((N-64)/138);
	sigma_ex=zeros(length(T1),1);
    for j=0:1:limit
        sigma_ex_flag=j*138+64;   
        for i=0.01:0.01:0.05
           x=i*v;
           sigma_ex(sigma_ex_flag,1)=(rp^2*acos((rp-x)/rp)-(rp-x)*(x*(2*rp-x))^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.06:0.01:0.1
           x=i*v;
           sigma_ex(sigma_ex_flag,1)=(rp^2*(pi-acos((x-rp)/rp))-(x-rp)*(x*(2*rp-x))^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.11:0.01:0.15
           x=i*v;
           sigma_ex(sigma_ex_flag,1)=(rp^2*(pi-acos((3*rp-x)/rp))-(3*rp-x)*(rp^2-(3*rp-x)^2)^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
        for i=0.16:0.01:0.2
           x=i*v;
           sigma_ex(sigma_ex_flag,1)=(rp^2*acos((x-3*rp)/rp)-(x-3*rp)*(rp^2-(x-3*rp)^2)^(1/2))*P/S0;
           sigma_ex_flag=sigma_ex_flag+1;
           if sigma_ex_flag>N
             break;
           end
        end
    end
end
%计算实际应力
function sigma=PlotRealSigma()
	global  w1 w2 T1 V1 alpha T
    star=1;
	N=length(T1);
	sigma=zeros(length(T1),1);
    for i=2:1:N
        temp=0;
        for j=star:1:i-1
            temp=temp+w2*(V1(j+1)+V1(j))*(T1(j+1)-T1(j))/2;
        end
        if mod(i,T)==0
            if star==1
                star=star+(T-1);
            else
                star=star+T;
            end
        end
        sigma(i,1)=temp+w1*V1(i);
    end
    sigma=-sigma./alpha;        
    
end
%读取文件
function [T1,V1]=ReadFile()
    
	% 2-2 滤波12.5 1M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1','B1:B1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1','C1:C1000');
	% 2-2 未滤波 1M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1','D1:D1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet1','E1:E1000');
	
	% 3-1 滤波12.5 1M
    %T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','B1:B1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','C1:C1000');
	% 3-1 未滤波 1M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','D1:D1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','E1:E1000');
	% 3-1 未滤波 2M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','G1:G1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','H1:H1000');
	% 3-1 滤波12.5 2M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','K1:K1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','L1:L1000');
	% 3-1 未滤波 4M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','L1:L1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','M1:M1000');
	% 3-1 未滤波 6M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','O1:O1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','P1:P1000');
	% 3-1 未滤波 8M
	T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','L1:L1000');
    V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet2','M1:M1000');
	
	% 4-2 未滤波 1M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','B1:B1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','C1:C1000');
	% 4-2 滤波12.5 1M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','D1:D1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','E1:E1000');
	% 4-2 未滤波 2M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','G1:G1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','H1:H1000');
    % 4-2 未滤波 4M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','L1:L1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','M1:M1000');
	% 4-2 未滤波 6M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','Q1:Q1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','R1:R1000');
	% 4-2 未滤波 8M
	%T1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','V1:V1000');
    %V1=xlsread('C:\\Users\\tx\\Desktop\\1.xlsx','Sheet3','W1:W1000');
end
