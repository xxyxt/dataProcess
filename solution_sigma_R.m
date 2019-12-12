clear;
clc;
syms xi z;
v=0.3;
h=0.065;
a=0.0178;
A=[1 1-2*v -1 1-2*v;
1 -2*v 1 2*v;
exp(-xi*h/a) (2-4*v+xi*h/a)*exp(-xi*h/a) exp(xi*h/a) -(2-4*v-xi*h/a)*exp(xi*h/a);
exp(-xi*h/a) -(2*v-xi*h/a)*exp(-xi*h/a) -exp(xi*h/a) (2*v+xi*h/a)*exp(xi*h/a)];
B=[1;0;0;0];
X=A\B;
% s(xi,z)=(X(1,1)+(1-2*v+xi*z/a)*X(2,1))*exp(-xi*z/a)-(X(3,1)-(1-2*v-xi*z/a)*X(4,1))*exp(xi*z/a);
fun_r(xi,z)=exp(-xi*z)*(X(1,1)-(1+2*v-xi*z)*X(2,1))-exp(xi*z)*(X(3,1)+(1+2*v+xi*z));
fun_sita(xi,z)=exp(-xi*z)*X(2,1)+exp(xi*z)*X(4,1);
fun_U=exp(-xi*z)*(X(1,1)-(1-xi*z)*X(2,1))-exp(xi*z)*(X(3,1)+(1+xi*z)*X(4,1));
N=500;
z_h=0.06;
step_z=0.001;
h=0.065;
step_x=0.1;
fun_array=zeros(N/step_x+1,z_h/step_z+1);
fun_array_r=zeros(N/step_x+1,z_h/step_z+1);
fun_array_sita=zeros(N/step_x+1,z_h/step_z+1);
fun_array_U=zeros(N/step_x+1,z_h/step_z+1);
index_y=1;
for z=0:step_z:z_h
    index_x=1;
    for xi=0:step_x:N
        temp_r=fun_r(xi,z);
        temp_sita=fun_sita(xi,z);
        temp_U=fun_U(xi,z);
        fun_array_r(index_x,index_y)=temp_r;
        fun_array_sita(index_x,index_y)=temp_sita;
        fun_array_U(index_x,index_y)=temp_U;
        index_x=index_x+1;
    end
    index_y=index_y+1;
end


