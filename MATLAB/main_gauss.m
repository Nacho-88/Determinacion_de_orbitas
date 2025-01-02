%DATOS CONCRETOS DEL PROBLEMA
deg=pi/180;         %deg2rad

mu = 398600;        %Parámetro gravitacional(kmˆ3/sˆ2)
RT = 6378;          %Radio Tierra (km)
f = 1/298.26;       %Factor de aplanamiento de la tierra
H =1;               %Elevación del puesto de observación (km)
phi = deg2rad(40);  %Latitud del puesto de observación (deg)

t=[0 118.104 237.577];   %Vector de tiempos de observación (seg)

AR = [43.5365 54.4196 64.3178]*deg;         %Vector de ascension recta (rad)
DEC = [-8.78334 -12.0739 -15.1054]*deg;     %Vector de declinación (rad)
theta= [44.5065 45.000 45.4992]*deg;        %Vector de LST (rad)

%Llamamos a la función que termina el problema

determinacion_orbita(t,AR,DEC,theta,mu,RT,phi,f,H)


