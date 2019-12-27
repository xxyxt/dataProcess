clear;
clc;
tic;
% [J_array,fun_array,res]=main();
% [limit,J_limit]=caculateLimit();
% [J_array,fun_array,res]=main_function();
% 计算函数矩阵
% fun_array=caculateFunArray();
% 计算积分
[res_array_r,res_array_sita]=cacalateInt();
toc;
function [res_array_r,res_array_sita]=cacalateInt()
    a=0.0178;
    path='D:\matlabworkbench\data';
    J_array_U=load('D:\matlabworkbench\data\J_array_U.mat');
    J_array_U=cell2mat(struct2cell(J_array_U));
    J_array=load('D:\matlabworkbench\data\J_array.mat');
    J_array=cell2mat(struct2cell(J_array));
    fun_array_U=load('D:\matlabworkbench\data\fun_array_U.mat');
    fun_array_U=cell2mat(struct2cell(fun_array_U));
    fun_array_sita=load('D:\matlabworkbench\data\fun_array_sita.mat');
    fun_array_sita=cell2mat(struct2cell(fun_array_sita));
    fun_array_r_new=load('D:\matlabworkbench\data\fun_array_r_new.mat');
    fun_array_r_new=cell2mat(struct2cell(fun_array_r_new));
    N=5000;
	res_array_U=zeros(60,1);
    res_array_r=zeros(60,1);
    res_array_sita=zeros(60,1);
    for j=1:1:60
		for i=1:1:N-1
			res_array_U(j)=res_array_U(j)+(fun_array_U(i+1,j)*J_array_U(i+1)+fun_array_U(i,j)*J_array_U(i))*a*0.1/2;
		end
    end
	for j=1:1:60
		for i=1:1:N-1
			res_array_r(j)=res_array_r(j)+(fun_array_r_new(i+1,j)*J_array(i+1)+fun_array_r_new(i,j)*J_array(i))*0.1/2;
        end
    end
	for j=1:1:60
		for i=1:1:N-1
			res_array_sita(j)=res_array_sita(j)+(fun_array_sita(i+1,j)*J_array(i+1)+fun_array_sita(i,j)*J_array(i))*2*0.3*0.1/2;
		end
    end
	res_array_r=res_array_U-res_array_r;
	res_array_sita=2*0.3*res_array_sita-res_array_U;
end
function fun_array=caculateFunArray()
    v=0.3;  % 泊松�?
    h=0.065;  % 板厚�?
    z=0;  % 深度
    a=0.0178;
    syms xi;
    fen_mu=v-2*exp(2*h*xi/a)-exp(4*h*xi/a)+8*v*exp(2*h*xi/a)+3*v*exp(4*h*xi/a)-8*v^2*exp(2*h*xi/a)+3*h*xi/a*exp(2*h*xi/a)+h*xi/a*exp(4*h*xi/a)-1;
    fen_zi_1=2*(v-v*exp(2*h*xi/a)-v^2+v^2*exp(2*h*xi/a)-h*xi/a*exp(2*h*xi/a)+v*h*xi/a*exp(2*h*xi/a));
    fen_zi_2=(2*v+xi/a*z-1)*(exp(2*h*xi/a)-v-v*exp(2*h*xi/a)+h*xi/a*exp(2*h*xi/a)+1);
    fen_zi_3=2*exp(2*h*xi/a)*(v+h*xi/a-v*exp(2*h*xi/a)-3*v^2+3*v^2*exp(2*h*xi/a)+v*xi/a*exp(2*h*xi/a));
    fen_zi_4=exp(2*h*xi/a)*(xi/a*z-2*v+1)*(v+exp(2*h*xi/a)-3*v*exp(2*h*xi/a)-h*xi/a*exp(2*h*xi/a)+1);
    f(xi)=exp(xi/a*z)*(fen_zi_1-fen_zi_2)/fen_mu-exp(-xi/a*z)*(fen_zi_3-fen_zi_4)/(fen_mu);
    N=500;
    step=0.1;
    index=1;
    fun_array=zeros(N/step+1,1);
    for xi=0:step:N
        fun_temp=f(xi);
        fun_array(index,1)=fun_temp;
        index=index+1;
    end
    % save('C:\Users\tx\Desktop\20191111\data\int\fun_array1.mat',fun_array);
end
function [J_array,fun_array,res]=main_function()
    v=0.3;  % 泊松�?
    h=0.065;  % 板厚�?
    z=0;  % 深度
    a=0.0178;
    syms xi;
    fen_mu=v-2*exp(2*h*xi/a)-exp(4*h*xi/a)+8*v*exp(2*h*xi/a)+3*v*exp(4*h*xi/a)-8*v^2*exp(2*h*xi/a)+3*h*xi/a*exp(2*h*xi/a)+h*xi/a*exp(4*h*xi/a)-1;
    fen_zi_1=2*(v-v*exp(2*h*xi/a)-v^2+v^2*exp(2*h*xi/a)-h*xi/a*exp(2*h*xi/a)+v*h*xi/a*exp(2*h*xi/a));
    fen_zi_2=(2*v+xi/a*z-1)*(exp(2*h*xi/a)-v-v*exp(2*h*xi/a)+h*xi/a*exp(2*h*xi/a)+1);
    fen_zi_3=2*exp(2*h*xi/a)*(v+h*xi/a-v*exp(2*h*xi/a)-3*v^2+3*v^2*exp(2*h*xi/a)+v*xi/a*exp(2*h*xi/a));
    fen_zi_4=exp(2*h*xi/a)*(xi/a*z-2*v+1)*(v+exp(2*h*xi/a)-3*v*exp(2*h*xi/a)-h*xi/a*exp(2*h*xi/a)+1);
    f(xi)=exp(xi/a*z)*(fen_zi_1-fen_zi_2)/fen_mu-exp(-xi/a*z)*(fen_zi_3-fen_zi_4)/(fen_mu);
    res=0;
    N=2000;
    step=0.1;
    index=1;
    J_array=zeros((N+1)/step,1);
    fun_array=zeros((N+1)/step,1);
    for xi=0:step:N
        A=besselj(1,xi+step);
        B=besselj(1,xi);
        fun_A=f(xi+step);
        fun_B=f(xi);
        res=res+(A*fun_A+B*fun_B)*step/2;
        J_array(index,1)=A;
        fun_array(index,1)=B;
        index=index+1;
    end
end
% 计算J积分的收敛点 返回�? limit=70
function [limit,J_limit]=caculateLimit()
    limit=0;
    alpha=1;
    J_limit=J(limit,alpha);
    while(abs(J(limit+1,alpha)-J_limit)>0.001)
        limit=limit+1;
        J_limit=J(limit,alpha);
    end

end
% 伽马函数计算 返回值gama_res
function gama_res=gama(m)
    syms t;
    gama_res=int(t^(m-1)/exp(t),t,0,inf);
end
% 计算J积分 输入参数 alpha为阶�? x为J积分的变�? 返回积分�?
function J_res=J(x,alpha)
        J_res=0;
        for m=0:1:500
            J_res_front=J_res;
            J_res=J_res+(-1)^m/(factorial(m)*gama(m+alpha+1))*(x/2)^(2*m+alpha);
            if abs(J_res_front-J_res)<0.000001  % 收敛�?
                break;
            end
        end
end
