function gibbs()

grados_rad=(2*pi)/360;

r=[50264855.39 52669175.89 55244867.63]; %Posición de mercurio en los dias 16,19 y 22 de diciembre 2024 (km)

long_helio=[133.9244 149.5687 163.7835]; %Longitud heliocéntrica de mercuio (grados)

lat_helio=[6.9833 6.8698 6.3279]; %Latitud heliocétrica de mercurio (grados)

long_helio_rad=long_helio.*grados_rad; %Lo pasamos a radianes

lat_helio_rad=lat_helio.*grados_rad;

mu=1.33e11; %Parámetro gravitacional Sol-Mercurio (km^3/s^2)

metodo_gibbs(r,long_helio_rad,lat_helio_rad,mu);

end