function gauss

grados_rad=(2*pi)/360;

mu=1.33e11;

r_T=6.957e5;

theta=[137.9711*grados_rad 138.3042*grados_rad 138.6368*grados_rad];

phi_prima=7.0034;

delta=[22.09813*grados_rad 21.99254*grados_rad 21.88647*grados_rad];

alpha=[142.72306*grados_rad 143.06115*grados_rad 143.39811*grados_rad];

t_1=0;

t_2=5400;

t_3=10800;

r_1x=r_T*cos(phi_prima)*cos(theta(1));
 
r_1y=r_T*cos(phi_prima)*sin(theta(1));
 
r_1z=r_T*sin(phi_prima);
 
r_1=[r_1x r_1y r_1z];
 
r_2x=r_T*cos(phi_prima)*cos(theta(2));
 
r_2y=r_T*cos(phi_prima)*sin(theta(2));

r_2z=r_T*sin(phi_prima);
 
r_2=[r_2x r_2y r_2z];
 
r_3x=r_T*cos(phi_prima)*cos(theta(3));
 
r_3y=r_T*cos(phi_prima)*sin(theta(3));

r_3z=r_T*sin(phi_prima);

r_3=[r_3x r_3y r_3z];

rho_1x=cos(delta(1))*cos(alpha(1));

rho_1y=cos(delta(1))*sin(alpha(1));

rho_1z=sin(delta(1));

rho_1=[rho_1x rho_1y rho_1z];

rho_2x=cos(delta(2))*cos(alpha(2));

rho_2y=cos(delta(2))*sin(alpha(2));

rho_2z=sin(delta(2));

rho_2=[rho_2x rho_2y rho_2z];

rho_3x=cos(delta(3))*cos(alpha(3));

rho_3y=cos(delta(3))*sin(alpha(3));

rho_3z=sin(delta(3));

rho_3=[rho_3x rho_3y rho_3z];

metodo_gauss(rho_1,rho_2,rho_3,r_1,r_2,r_3,t_1,t_2,t_3,mu)

end