function metodo_gibbs(r,long_helio,lat_helio,mu)

r_1=r(1); %Observación 1 

r_2=r(2); %Observación 2

r_3=r(3); %Observación 3

lh_1=long_helio(1);

lh_2=long_helio(2);

lh_3=long_helio(3);

lah_1=lat_helio(1);

lah_2=lat_helio(2);

lah_3=lat_helio(3);

%Calculamos el vector posición.
r1_vect=[r_1*cos(lah_1)*cos(lh_1) r_1*cos(lah_1)*sin(lh_1) r_1*sin(lah_1)];

r2_vect=[r_2*cos(lah_2)*cos(lh_2) r_2*cos(lah_2)*sin(lh_2) r_2*sin(lah_2)];

r3_vect=[r_3*cos(lah_3)*cos(lh_3) r_1*cos(lah_3)*sin(lh_3) r_1*sin(lah_3)];

c23=cross(r2_vect,r3_vect);

a=dot(r1_vect,c23)/(r_1*norm(c23));

S=r1_vect.*(r_2-r_3)+r2_vect.*(r_3-r_1)+r3_vect.*(r_1-r_2);

D=cross(r1_vect,r2_vect)+cross(r2_vect,r3_vect)+cross(r3_vect,r1_vect);

N=r_1*cross(r2_vect,r3_vect)+r_2*cross(r3_vect,r1_vect)+r_3*cross(r1_vect,r2_vect);

v=zeros(3,3); %Vector de velocidades (km/s)

v(1,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r1_vect)/r_1+S);

v(2,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r2_vect)/r_2+S);

v(3,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r3_vect)/r_3+S);

v_1=v(1,:); %Vectror velocidad para la observación 1 

v_2=v(2,:); %Vector de velocidad para la observación 2 

v_3=v(3,:); %Vector de velocidad para la observación 3

fprintf('Seleccione qué vectores de posición y velocidad quiere usar para calcular los parámetros orbitales:\n')
var=input('1 para observación 1. 2 para observación 2 y 3 para observación 3\n');

if var==1 || var==2 || var==3

    if var==1

    calc_elem_orb(r1_vect,v_1,mu)

    end

    if var==2 

        calc_elem_orb(r2_vect,v_2,mu)

    end

    if var==3 

        calc_elem_orb(r3_vect,v_3,mu)

    end

else

    fprintf('No ha seleccionado un número válido')
    return;
    
end