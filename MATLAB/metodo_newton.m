function [sol,raiz] = metodo_newton(f,x0,epsilon,epsilon_2,df)

N=100;

n=0;

x=zeros(1,N+1);

x(1)=x0;

for k=1: N

    x(k+1)=x(k)-(f(x(k))/df(x(k)));

    n=n+1;

    if (abs(x(k+1)-x(k))/x(k)<epsilon)&&(abs(f(x(k+1)))<epsilon_2)

        x(k+2:end)=[];

        sol=x;

        raiz=x(k+1);

        fprintf("El metodo ha convergido con %d iteraciones\n",n)

        return

    end

end

fprintf("Se ha alcanzado el numero maximo de iteraciones y aun no converge\n")

raiz=x(end);

sol=x;

end