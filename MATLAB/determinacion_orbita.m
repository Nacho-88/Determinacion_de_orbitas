function determinacion_orbita(t,AR,DEC,theta,mu,RT,phi,f,H)

fac1 = RT/sqrt(1-(2*f- f*f)*sin(phi)^2);
fac2 = (RT*(1-f)^2/sqrt(1-(2*f- f*f)*sin(phi)^2) + H) *sin(phi);

for i = 1:3
 R(i,1) = (fac1 + H)*cos(phi)*cos(theta(i));
 R(i,2) = (fac1 + H)*cos(phi)*sin(theta(i));
 R(i,3) = fac2;
 rho(i,1) = cos(DEC(i))*cos(AR(i));
 rho(i,2) = cos(DEC(i))*sin(AR(i));
 rho(i,3) = sin(DEC(i));
 end

[r, v, r_old, v_old] = metodo_gauss(rho(1,:), rho(2,:), rho(3,:), R(1,:), ...
R(2,:), R(3,:), t(1), t(2), t(3),mu);

fprintf('Resultados sin mejora:\n')
fprintf('\nr (km)= [%g, %g, %g]',r_old(1), r_old(2), r_old(3))
fprintf('\nv (km/s)= [%g, %g, %g]',v_old(1), v_old(2), v_old(3))
calc_elem_orb(r_old,v_old,mu);

fprintf('Resultados con la mejora iterativa:\n')
fprintf('\nr (km)= [%g, %g, %g]',r(1), r(2), r(3))
fprintf('\nv (km/s)= [%g, %g, %g]',v(1), v(2), v(3))
calc_elem_orb(r,v,mu);

end
