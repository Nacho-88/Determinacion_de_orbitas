function metodo_gauss(rho_1,rho_2,rho_3,r_1,r_2,r_3,t_1,t_2,t_3,mu)

tau_1=t_1-t_2;

tau_3=t_3-t_2;

tau=tau_3-tau_1;

p_1=cross(rho_2,rho_3);

p_2=cross(rho_1,rho_3);

p_3=cross(rho_1,rho_2);

D_0=dot(rho_1,p_1);

D_11=dot(r_1,p_1);

D_21=dot(r_2,p_1);

D_31=dot(r_3,p_1);

D_12=dot(r_1,p_2);

D_22=dot(r_2,p_2);

D_32=dot(r_3,p_2);

D_13=dot(r_1,p_3);

D_23=dot(r_2,p_3);

D_33=dot(r_3,p_3);

A=(1/D_0)*(-D_12*(tau_3/tau)+D_22+D_32*(tau_1/tau));

B=1/(6*D_0)*(D_12*(tau_3^2-tau^2)*(tau_3/tau)+D_32*(tau^2-tau_1^2)*(tau_1/tau));

E=dot(r_2,rho_2);

r_2sq=dot(r_2,r_2);

a=-(A^2+2*A*E+r_2sq);

b=-2*mu*B*(A+E);

c=-mu^2*B^2;

Roots = roots([1 0 a 0 0 b 0 0 c]);

x=Roots(Roots>0 & ~imag(Roots) );

aux=length(x);

if aux>1

    times=[-1e4,100,1e4];

    func= @(x) x.^8+a*x.^6+b*x.^3+c;

    plot(times,func(times))

    xlabel('x')

    ylabel('F')

   fprintf('Tienes que quedarte con solo 1 raíz. Hay %d raíces\n',aux)

   fprintf('Esas raíces son [%s]\n', sprintf('%.2f, ', x(1:end))')

   aux_2=input('Selecciona con cuál quieres quedarte para el método con los números de tu teclado\n');

   x=x(aux_2);

   if aux_2 < 1 || aux_2>aux

       fprintf('No has seleccionado un número válido. El programa va a finalizar\n')

       return;

   end

end

ayuda=6*(D_31*(tau_1/tau_3)+D_21*(tau/tau_3))*x^3+mu*D_31*(tau^2-tau_1^2)*(tau_1/tau_3);

ayuda_2=6*x^3+mu*(tau^2-tau_3^2);

Rho_1=1/(D_0)*(ayuda/(ayuda_2)-D_11);

ayuda_3=6*(D_13*(tau_3/tau_1)-D_23*(tau/tau_1))*x^3+mu*D_13*(tau^2-tau_3^2)*(tau_3/tau_1);

ayuda_4=6*x^3+mu*(tau^2-tau_3^2);

Rho_3=1/(D_0)*(ayuda_3/(ayuda_4)-D_33);

Rho_2=A+(mu*B)/x^3;

R_1=r_1+Rho_1.*rho_1;

R_2=r_2+Rho_2.*rho_2;

R_3=r_3+Rho_3.*rho_3;

f_1=1-0.5*(mu/x^3)*tau_1^2;

f_3=1-0.5*(mu/x^3)*tau_3^2;

g_1=tau_1-(1/6)*(mu/x^3)*tau_1^3;

g_3=tau_3-(1/6)*(mu/x^3)*tau_3^3;

auux=1/(f_1*g_3-f_3*g_1);

auux_2=-f_3.*R_1+f_1.*R_3;

v_2=auux.*auux_2;

calc_elem_orb(R_2,v_2,mu)

end