clear;
clc;
tic;
syms xi z;
v=0.3;
h=0.05;
a=0.0178;
% 解方程
A=[1 1-2*v -1 1-2*v;
1 -2*v 1 2*v;
exp(-xi*h/a) (2-4*v+xi*h/a)*exp(-xi*h/a) exp(xi*h/a) -(2-4*v-xi*h/a)*exp(xi*h/a);
exp(-xi*h/a) -(2*v-xi*h/a)*exp(-xi*h/a) -exp(xi*h/a) (2*v+xi*h/a)*exp(xi*h/a)];
B=[1;0;0;0];
X=A\B;
% 函数表达式
fun_z(xi,z)=(X(1,1)+(1-2*v+xi*z/a)*X(2,1))*exp(-xi*z/a)-(X(3,1)-(1-2*v-xi*z/a)*X(4,1))*exp(xi*z/a);
fun_r(xi,z)=exp(-xi*z/a)*(X(1,1)-(1+2*v-xi*z/a)*X(2,1))-exp(xi*z/a)*(X(3,1)+(1+2*v+xi*z/a)*X(4,1));
fun_sita(xi,z)=exp(-xi*z/a)*X(2,1)+exp(xi*z/a)*X(4,1);
fun_U(xi,z)=exp(-xi*z/a)*(X(1,1)-(1-xi*z/a)*X(2,1))-exp(xi*z/a)*(X(3,1)+(1+xi*z/a)*X(4,1));
% 初始化函数矩阵
N=2000;
z_h=0.05;
step_z=0.001;
step_x=0.1;
fun_array_z=zeros(N/step_x+1,z_h/step_z+1);
fun_array_r=zeros(N/step_x+1,z_h/step_z+1);
fun_array_sita=zeros(N/step_x+1,z_h/step_z+1);
fun_array_U=zeros(N/step_x+1,z_h/step_z+1);
 % 列下标
index_y=1;
% 代入不同深度求函数值 这个循环比较耗时间 
for z=0:step_z:z_h
    % 行下标
    index_x=1; 
    for xi=0:step_x:N
        temp_z=fun_z(xi,z);
        % temp_r=fun_r(xi,z);
        % temp_sita=fun_sita(xi,z);
        % temp_U=fun_U(xi,z);

        fun_array_z(index_x,index_y)=temp_z;
        % fun_array_r(index_x,index_y)=temp_r;
        % fun_array_sita(index_x,index_y)=temp_sita;
        % fun_array_U(index_x,index_y)=temp_U;
        index_x=index_x+1;
    end
    index_y=index_y+1;
end
% 求J积分
J_array=zeros(N/step_x+1,1);
index_J=1;
for xi=0:0.1:N
    J_array(index_J,1)=besselj(1,xi);
    index_J=index_J+1;
end
% 初始化res矩阵
res_array_z=zeros(z_h/step_z,1);
res_array_r=zeros(z_h/step_z,1);
res_array_sita=zeros(z_h/step_z,1);
res_array_U=zeros(z_h/step_z,1);
% 数值积分
for j=1:1:z_h/step_z
    for i=1:1:N/step_x
        res_array_z(j,1)=res_array_z(j,1)+(fun_array_z(i+1,j)*J_array(i+1)+fun_array_z(i,j)*J_array(i))*0.1/2;
        % res_array_r(j,1)=res_array_r(j,1)+(fun_array_r(i+1,j)*J_array(i+1)+fun_array_r(i,j)*J_array(i))*0.1/2;
        % res_array_sita(j,1)=res_array_sita(j,1)+2*v*(fun_array_sita(i+1,j)*J_array(i+1)+fun_array_sita(i,j)*J_array(i))*0.1/2;
        % res_array_U(j,1)=res_array_U(j,1)+0.5*(fun_array_U(i+1,j)*J_array(i+1)+fun_array_U(i,j)*J_array(i))*0.1/2;
    end
end
res_array_r=res_array_U-res_array_r;
res_array_sita=0.6*res_array_sita-res_array_U;
toc;
% 写成函数形式
% global fun_array_z fun_array_r fun_array_sita fun_array_U N z_h step_x step_z v h a X
% v=0.3;
% h=0.05;
% a=0.0178;
% N=1000;
% z_h=0.05;
% step_z=0.001;
% step_x=0.1;
% fun_array_z=zeros(N/step_x+1,z_h/step_z+1);
% fun_array_r=zeros(N/step_x+1,z_h/step_z+1);
% fun_array_sita=zeros(N/step_x+1,z_h/step_z+1);
% fun_array_U=zeros(N/step_x+1,z_h/step_z+1);
% X=dealSolution();
% % loop();
% function loop()  
%     global N z_h step_z
%     % 代入不同深度求函数值
%     for z=0:step_z:z_h
%         index_y=z/0.001+1;
%         parfor xi=0:1:N/0.1
%             index_x=xi+1;
            
%             loopBlock(xi,z,index_x,index_y);
%         end
%     end
% end
% function loopBlock(xi,z,index_x,index_y)
%     global  fun_array_z fun_array_r fun_array_sita fun_array_U step_x
%     % 初始化下标
%     % index_y=z/step_z+1;
%     % index_x=xi+1;
%     % 代入计算
%     temp_z=fun_z(xi*step_x,z);
%     temp_r=fun_r(xi*step_x,z);
%     temp_sita=fun_sita(xi*step_x,z);
%     temp_U=fun_U(xi*step_x,z);
%     % 赋值
%     fun_array_z(index_x,index_y)=temp_z;
%     fun_array_r(index_x,index_y)=temp_r;
%     fun_array_sita(index_x,index_y)=temp_sita;
%     fun_array_U(index_x,index_y)=temp_U;
% end
% function X=dealSolution()
%     global X v h a
%     syms xi z;
%     % 解方程
%     A=[1 1-2*v -1 1-2*v;
%     1 -2*v 1 2*v;
%     exp(-xi*h/a) (2-4*v+xi*h/a)*exp(-xi*h/a) exp(xi*h/a) -(2-4*v-xi*h/a)*exp(xi*h/a);
%     exp(-xi*h/a) -(2*v-xi*h/a)*exp(-xi*h/a) -exp(xi*h/a) (2*v+xi*h/a)*exp(xi*h/a)];
%     B=[1;0;0;0];
%     X=A\B;
% end
% function res_fun_z=fun_z(xi,z)
%     global X v a
%     res_fun_z=(X(1,1)+(1-2*v+xi*z/a)*X(2,1))*exp(-xi*z/a)-(X(3,1)-(1-2*v-xi*z/a)*X(4,1))*exp(xi*z/a);
% end
% function res_fun_r=fun_r(xi,z)
%     global X v a
%     res_fun_r=exp(-xi*z/a)*(X(1,1)-(1+2*v-xi*z/a)*X(2,1))-exp(xi*z/a)*(X(3,1)+(1+2*v+xi*z/a)*X(4,1));
% end
% function res_fun_sita=fun_sita(xi,z)
%     global X a
%     res_fun_sita=exp(-xi*z/a)*X(2,1)+exp(xi*z/a)*X(4,1);
% end
% function res_fun_U=fun_U(xi,z)
%     global X a
%     res_fun_U=exp(-xi*z/a)*(X(1,1)-(1-xi*z/a)*X(2,1))-exp(xi*z/a)*(X(3,1)+(1+xi*z/a)*X(4,1));
% end
