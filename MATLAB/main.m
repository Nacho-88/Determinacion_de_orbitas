fprintf(['En este script puedes seleccionar que método usar para calcular ' ...
    'la órbita de Mercurio.\nAdicionalmente, también puedes usar el ' ...
    'método de Gauss para calcular la órbita del satélite Meteosat 10\n']);
a=input(['Para calcular la órbita de Mercurio por el método de Gibbs pulse ' ...
    '1.\nPara calcular la órbita de Mercurio por el método de Gauss pulse ' ...
    '2.\nPara calcular la órbita del Meteosat 10 pulse 3\n']);

if a==1

    gibbs();

end

if a==2

    gauss();

end

if a==3

    gauss_sat(); 

end

if a~=1 && a~=2 && a~=3 

    fprintf(['No ha introducido un número correcto. Por favor seleccione 1 ' ...
        '2 y 3.\nTiene que volver a correr el programa para ello.\n'])

end