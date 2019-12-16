syms q E3
E=6.5*10^10;
v=0.3;
F_r=-0.24336*q;
F_z=-0.93164*q;
c11=E;
c12=-v*E;
c13=-v*E;
c33=E;
d31=275*10^-12;
d33=650*10^-12;
epslon_x=(1/2^0.5)*(F_r-v*F_z)/E;
epslon_y=(1/2^0.5)*(F_r-v*F_z)/E;
epslon_z=(F_z-v*F_r)/E;
mu_33=3850*8.854*10^-12;
sigma_x=(c11*epslon_x+c12*epslon_y+c13*epslon_z)-d31*E3;
sigma_y=(c12*epslon_x+c11*epslon_y+c13*epslon_z)-d31*E3;
sigma_z=(c13*epslon_x+c13*epslon_y+c33*epslon_z)-d31*E3;
D_z=d31*sigma_x+d31*sigma_y+d33*sigma_z+mu_33*E3;

