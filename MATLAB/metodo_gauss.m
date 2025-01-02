%-------------------------------------------------------------------------
function [r, v, r_old, v_old] = metodo_gauss(Rho1, Rho2, Rho3, R1, R2, R3, t1, t2, t3,mu)

% This function uses the Gauss method with iterative
% improvement (Algorithms 5.5 and 5.6) to calculate the state
% vector of an orbiting body from angles-only observations at
% three closely-spaced times.

% mu                - the gravitational parameter (km^3/s^2)

% t1, t2, t3        - the times of the observations (s)

% tau, tau1, tau3   - time intervals between observations (s)

% R1, R2, R3        - the observation site position vectors
%                     at t1, t2, t3 (km)

% Rho1, Rho2, Rho3  - the direction cosine vectors of the
%                     satellite at t1, t2, t3

% p1, p2, p3        - cross products among the three direction
%                     cosine vectors

% Do                - scalar triple product of Rho1, Rho2 an Rho3

% D                 - Matrix of the nine scalar triple products
%                     of R1, R2 and R3 with p1, p2 and p3

% E                 - dot product of R2 and Rho2

% A, B              - constants in the expression relating
%                     slant range to geocentric radius

% a,b,c             - coefficients of the 8th order polynomial
%                     in the estimated geocentric radius x

% x                 - positive root of the 8th order polynomial

% rho1, rho2, rho3  - the slant ranges at t1, t2, t3

% r1, r2, r3        - the position vectors at t1, t2, t3 (km)

% r_old, v_old      - the estimated state vector at the end of
%                     Algorithm 5.5 (km, km/s)

% rho1_old,
% rho2_old, and
% rho3_old          - the values of the slant ranges at t1, t2,
%                     t3 at the beginning of iterative
%                     improvement (Algorithm 5.6) (km)

% diff1, diff2,
% and diff3         - the magnitudes of the differences between
%                     the old and new slant ranges at the end
%                     of each iteration

% tol               - the error tolerance determining convergence
%

% n                 - number of passes through the iterative improvement loop

% nmax              - limit on the number of iterations

% ro, vo            - magnitude of the position and velocity vectors (km, km/s)

% vro               - radial velocity component (km)

% alpha             - reciprocal of the semimajor axis (1/km)

% v2                - computed velocity at time t2 (km/s)

% r, v              - the state vector at the end of Algorithm 5.6 (km, km/s)

%-------------------------------------------------------------------------

%Ecuaciones 4.7 , 4.8 y 4.9:
tau1 = t1- t2;
tau3 = t3- t2;

tau = tau3- tau1;

%Ecuaciones 4.10:

p1 = cross(Rho2,Rho3);
p2 = cross(Rho1,Rho3);
p3 = cross(Rho1,Rho2);

%Ecuación 4.11
Do = dot(Rho1,p1);

%Ecuaciones 4.12
D = [[dot(R1,p1) dot(R1,p2) dot(R1,p3)]
    [dot(R2,p1) dot(R2,p2) dot(R2,p3)]
    [dot(R3,p1) dot(R3,p2) dot(R3,p3)]];

%Cálculo de E
E = dot(R2,Rho2);

%Ecuaciones 4.26 y 4.27
A = 1/Do*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);

B = 1/6/Do*(D(1,2)*(tau3^2- tau^2)*tau3/tau + D(3,2)*(tau^2- tau1^2)*tau1/tau);

%Coeficientes polinomio grado 8
a =-(A^2 + 2*A*E + norm(R2)^2);
b =-2*mu*B*(A + E);
c =-(mu*B)^2;

%Calcular las raices del polinomio de la ecuación 4.32 usando la fc de MATLAB
Roots = roots([1 0 a 0 0 b 0 0 c]);

%Encontramos la raiz positiva real
raiz = posroot(Roots);

%Ecuaciones 4.17
f1 = 1- 1/2*mu*tau1^2/raiz^3;
f3 = 1- 1/2*mu*tau3^2/raiz^3;

g1 = tau1- 1/6*mu*(tau1/raiz)^3;
g3 = tau3- 1/6*mu*(tau3/raiz)^3;

%Ecuación 4.25
rho2=A+mu*B/raiz^3;

%Ecuaciones 4.28 y 4.29
rho1 =  1/Do*((6*(D(3,1)*tau1/tau3 + D(2,1)*tau/tau3)*raiz^3 ...
        + mu*D(3,1)*(tau^2- tau1^2)*tau1/tau3) ...
        /(6*raiz^3 + mu*(tau^2- tau3^2))- D(1,1));

rho3 =  1/Do*((6*(D(1,3)*tau3/tau1- D(2,3)*tau/tau1)*raiz^3 ...
        + mu*D(1,3)*(tau^2- tau3^2)*tau3/tau1) ...
        /(6*raiz^3 + mu*(tau^2- tau3^2))- D(3,3));

%Ecuaciones 4.1, 4.2 y 4.3
r1 = R1 + rho1*Rho1;
r2 = R2 + rho2*Rho2;
r3 = R3 + rho3*Rho3;

%Ecuación 4.35
v2 = (-f3*r1 + f1*r3)/(f1*g3- f3*g1);

%Guardamos las primeras estimaciones de r2 y v2
r_old = r2;
v_old = v2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AQUI ACABARÍA EL PRIMER ALGORITMO DE DETERMINACIÓN, AHORA PASAMOS A ITERAR
%LO QUE TENEMOS PARA AUMENTAR LA PRECISIÓN DE LOS RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inicializamos el bucle iterativo de mejora y establecemos la tolerancia
rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;

diff1=1;    diff2=1;    diff3=1;

n =0;
nmax = 1000;
tol = 1e-8;

%Inicio del bucle
while ((diff1 > tol) && (diff2 > tol) && (diff3 > tol)) && (n < nmax)

    n = n+1;

    %Datos necesarios para las ecuaciones de Kepler
    ro = norm(r2);
    vo = norm(v2);
    vro = dot(v2,r2)/ro;
    alpha = 2/ro- vo^2/mu;

    %Resolvemos  la ecuacion universal de kepler para los tiempos t1 y t3
    x1 = kepler_U(tau1, ro, vro, alpha,mu);
    x3 = kepler_U(tau3, ro, vro, alpha,mu);

    %Calculamos los coeficientes de Lagrange para los tiempos t1 y t3
    [ff1, gg1] = f_and_g(x1, tau1, ro, alpha,mu);
    [ff3, gg3] = f_and_g(x3, tau3, ro, alpha,mu);

    %Ahora vamos a actualizar los valores de f y g haciendo un promedio entre
    %los que teníamos y los que tenemos ahora.
    f1 = (f1 + ff1)/2;
    f3 = (f3 + ff3)/2;
    g1 = (g1 + gg1)/2;
    g3 = (g3 + gg3)/2;

    %Ecuaciones 4.15 y 4.16
    c1 = g3/(f1*g3- f3*g1);
    c3 =-g1/(f1*g3- f3*g1);

    %Ecuaciones 4.22, 4.23 y 4.24
    rho1 = 1/Do*(-D(1,1) + 1/c1*D(2,1)- c3/c1*D(3,1));
    rho2 = 1/Do*(-c1*D(1,2) + D(2,2)- c3*D(3,2));
    rho3 = 1/Do*(-c1/c3*D(1,3) + 1/c3*D(2,3)- D(3,3));

    %Y ahora usamos otra vez las ecuaciones 4.1, 4.2 y 4.3 para actualizar
    %los valores de las posiciones ahora que ya las hemos mejorado
    r1 = R1 + rho1*Rho1;
    r2 = R2 + rho2*Rho2;
    r3 = R3 + rho3*Rho3;

    %Así como volvemos a usar la ecuación 4.35
    v2 = (-f3*r1 + f1*r3)/(f1*g3- f3*g1);

    %Acabando, calculamos las diferencias en los valores rho para ver como 
    %va el criterio de convergencia
    diff1 = abs(rho1- rho1_old);
    diff2 = abs(rho2- rho2_old);
    diff3 = abs(rho3- rho3_old);

    %y por último establecemos como antiguos los que acabamos de calcular
    rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;

end
%Fin del bucle

fprintf('\n( **Numero de iteraciones del bucle de mejora')
fprintf(' = %g)\n\n', n)

if n >= nmax
    fprintf('\n\n **El número de iteraciones del bucle supera el límite de %g \n\n ', nmax);
end

%Por último devolvemos el vector de estados para la observación intermedia

 r = r2;
 v = v2;



end
