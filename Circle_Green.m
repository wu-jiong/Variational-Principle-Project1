clear;clc;
x=linspace(-1,1,100);
y=linspace(-1,1,100);
dx=x(2)-x(1);
dy=y(2)-y(1);
Nx=length(x);
Ny=length(y);

f0=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        f0(i,j)=(x(i)).^2+(y(j)).^2;
    end
end
sigma=0.1;
f0=-f0/(2*sigma^2);
f0=exp(f0)/(2*pi*sigma^2);
% f0(50:70,30:50)=1;  % Another initial trial function
surf(x,y,f0);
title('Iterations k=0',FontSize=14)

G= @(x,y,xp,yp) 1/4/pi*log((x^2+xp^2-2*x*xp*cos(y-yp))/(x^2*xp^2+1-2*x*xp*cos(y-yp)));

lambda=1;
Ls=zeros(10,1);
Es=zeros(10,1);

for k=1:10
    f0=lambda*f0;
    f1=zeros(Nx,Ny);
    for i=1:Nx
        for j=1:Ny
            rho=sqrt((x(i))^2+(y(j))^2);
            if rho>1
                continue
            end

            if y(j)>=0
                theta=acos(x(i)/rho);
            end
            if y(j)<0
                theta=2*pi-acos(x(i)/rho);
            end

            for m=1:Nx
                for n=1:Ny
                    rhop=sqrt((x(m))^2+(y(n))^2);
                    if rhop>1
                        continue
                    end

                    if y(n)>=0
                        thetap=acos(x(m)/rhop);
                    end
                    if y(n)<0
                        thetap=2*pi-acos(x(m)/rhop);
                    end

                    if rhop==rho&&thetap==theta
                        continue
                    end
                    f1(i,j)=f1(i,j)-dx*dy*G(rho,theta,rhop,thetap)*(f0(m,n)-tanh(f0(m,n)).^2);
                end
            end

        end
    end
    f0=f1;

    % Solve Lambda
    dfx=zeros(Nx,Ny);
    dfy=zeros(Nx,Ny);
    for i=1:Nx-1
        dfx(i,:)=(f0(i+1,:)-f0(i,:))/dx;
    end
    dfx(Nx,:)=dfx(Nx-1,:);
    for j=1:Ny-1
        dfy(:,j)=(f0(:,j+1)-f0(:,j))/dy;
    end
    dfy(:,Ny)=dfy(:,Ny-1);
    
    num=0;
    den=0;
    for i=1:Nx
        for j=1:Ny
            r_g=sqrt(x(i)^2+y(j)^2);
            if r_g>1
                continue
            end
            den=den+f0(i,j)^3*dx*dy;
            num=num+(f0(i,j)^2+dfx(i,j)^2+dfy(i,j)^2)*dx*dy;
        end
    end
    lambda=num/den;
    Ls(k)=lambda;

    error=check_error(f0,x,y,1);
    a=0;
    E=0;
    for i=1:Nx
        for j=1:Nx
            if x(i)^2+y(j)^2<0.9
                a=a+1;
                E=E+error(i,j);
            end
        end
    end
    Es(k)=E/a;

    figure(2)
    surf(x,y,f0)
    title(['Interations k=',num2str(k)],FontSize=14)
    
end




