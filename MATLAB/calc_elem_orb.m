function calc_elem_orb(R,V,mu)

rad_grados=360/(2*pi);

r=norm(R);

v=norm(V);

v_r=dot(R,V)/r;

H=cross(R,V);

h=norm(H);

fprintf('\nEl valor del momento angular específico es de %.2e km^2/s\n',h)

i=acos(H(3)/h);

fprintf('El valor de la inclinación es de %f º\n',i*rad_grados)

N = cross([0 0 1],H);

n = norm(N);

if n ~= 0

R_A = acos(N(1)/n);

if N(2) < 0

R_A = 2*pi - R_A;

end

else

R_A = 0;

end

fprintf('El valor del nodo ascendente es de %f º\n',R_A*rad_grados)

E = 1/mu*((v^2 - mu/r).*R - r*v_r.*V);

e = norm(E);

fprintf('El valor de la excentricidad es de %f\n',e)

w = acos(dot(N,E)/n/e);

if E(3) < 0

    w = 2*pi - w;

end

fprintf('El valor del argumento del perigeo es de %f º\n',w*rad_grados)

a = h^2/(mu*(1 - e^2));

fprintf('El valor del semieje mayor es de %.2e km\n',a)

theta= acos(dot(E,R)/(e*r));

if v_r < 0

    theta= 2*pi - theta;

end

fprintf('El valor de la anomalía verdadera es de %f º\n',theta*rad_grados)

fprintf('Radio del periastro = %g km\n', h^2/mu/(1 + e))


if e < 1

    T = 2*pi*sqrt(a^3/mu);
    fprintf('Periodo:')
    fprintf('\nSegundos = %g', T)
    fprintf('\nMinutos = %g', T/60)
    fprintf('\nHoras = %g', T/3600)
    fprintf('\nDías = %g\n\n', T/24/3600)


end


end





