function p = metodo_gibbs(T,e,a,mu,t_trans)

UA_m=149597870700;

horas_seg=86400;

p=zeros(1,3); %Semilatus rectum

p(1)=a(1)*(1-e(1)^2); 

p(2)=a(2)*(1-e(2)^2);

p(3)=a(3)*(1-e(3)^2);

theta=zeros(1,3); %Anomalía verdadera

M=zeros(1,3);

E0=0;

epsilon=1e-12;

epsilon_2=epsilon;

for i=1 : 3

    t_T0=t_trans*horas_seg;

    n=(2*pi)/T(i);

    M(i)=n*t_T0;

    func= @(E) M(i)-E+e(i)*sin(E);

    d_func= @(E) -1+e(i)*cos(E);

    [~,E]=metodo_newton(func,E0,epsilon,epsilon_2,d_func);

    theta(i)=2*atan(sqrt((1+e(i))/(1-e(i)))*tan(E/2));

end

r=zeros(1,3); %Vector de posición (UA)

r(1)=p(1)/(1+e(1)*cos(theta(1)));

r(2)=p(2)/(1+e(2)*cos(theta(2)));

r(3)=p(3)/(1+e(3)*cos(theta(3)));

r_m=zeros(1,3); %Vector posicion en metros

r_m(1)=r(1)*UA_m;

r_m(2)=r(2)*UA_m;

r_m(3)=r(3)*UA_m;

r_Sp=zeros(3,3); %Vector posición en el sistema perifocal

r_Sp(1,:)= [r_m(1)*cos(theta(1)), r_m(1)*sin(theta(1)),0];

r_Sp(2,:)= [r_m(2)*cos(theta(2)), r_m(2)*sin(theta(2)),0];

r_Sp(3,:)= [r_m(3)*cos(theta(3)), r_m(3)*sin(theta(3)),0];

S=r_Sp(1,:).*(r(2)-r(3))+r_Sp(2,:).*(r(3)-r(1))+r_Sp(3,:).*(r(1)-r(2));

D=cross(r_Sp(1,:),r_Sp(2,:))+cross(r_Sp(2,:),r_Sp(3,:))+cross(r_Sp(3,:),r_Sp(1,:));

N=norm(r_Sp(1,:))*cross(r_Sp(2,:),r_Sp(3,:))+norm(r_Sp(2,:))*cross(r_Sp(3,:),r_Sp(1,:))+norm(r_Sp(3,:))*cross(r_Sp(1,:),r_Sp(2,:));

v=zeros(3,3); %Vector de velocidades

v(1,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r_Sp(1,:))/r_m(1)+S);

v(2,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r_Sp(2,:))/r_m(2)+S);

v(3,:)=sqrt(mu/(norm(N)*norm(D)))*(cross(D,r_Sp(3,:))/r_m(3)+S);

p=[r_Sp; v];

end