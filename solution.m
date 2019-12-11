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
s(xi,z)=(X(1,1)+(1-2*v+xi*z/a)*X(2,1))*exp(-xi*z/a)-(X(3,1)-(1-2*v-xi*z/a)*X(4,1))*exp(xi*z/a);
N=500;
z_h=0.06;
step_z=0.001;
h=0.065;
step_x=0.1;
fun_array=zeros(N/step_x+1,z_h/step_z+1);
index_y=1;
for z=0.041:step_z:z_h
    index_x=1;
    for xi=0:step_x:N
        temp=s(xi,z);
        fun_array(index_x,index_y)=temp;
        index_x=index_x+1;
    end
    index_y=index_y+1;
end
path='C:\Users\tx\Desktop\20191111\data\int';
fileName='fun_array_different_z.xlsx';
xlswrite([path fileName],fun_array);


