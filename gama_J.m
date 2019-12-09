i=1;
J_array=zeros(51,1);
x_array=zeros(51,1);
for x=0:1:50
    J_res=J(x,0);
    J_array(i,1)=J_res;
    x_array(i,1)=x;
    i=i+1;
end
plot(x_array,J_array);
function gama_res=gama(m)
    syms t;
    gama_res=int(t^(m-1)/exp(t),t,0,inf);
end
function J_res=J(x,alpha)
    J_res=0;
    for m=0:1:500
        J_res=J_res+(-1)^m/(factorial(m)*gama(m+alpha+1))*(x/2)^(2*m+alpha);
    end
end
