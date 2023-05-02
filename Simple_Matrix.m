clear;clc

h=0.02;    
a=1;     
b=1;     
m=2*a/h;   
n=2*b/h;     

x = -a:h:a;
y = -b:h:b;
Nx=length(x);
Ny=length(y);
dx=h;
dy=h;
f=zeros(Nx,Ny); 

for i=1:Nx
    for j=1:Ny
        f(i,j)=(x(i)).^2+(y(j)).^2;
    end
end
sigma=0.1;
f=-f/(2*sigma^2);
f=exp(f)/(2*pi*sigma^2);
surf(x,y,f)
lambda=1;

q=1.5;
for kk=1:10
    s=1; 
    sum=0; 
    v=zeros(n+1,m+1);
    while s>1e-6   
        w=v;   
        for k=2:1:n
            for l=2:1:m
                rho=sqrt((x(k))^2+(y(l))^2);
                if rho>1
                    continue
                end
                v(k,l)=(1-q)*v(k,l)+q*(v(k-1,l)+v(k,l-1)+v(k+1,l)+v(k,l+1)-h*h*(f(k,l)^2))/4;
            end
        end
        s=max(max(abs(v-w)));    
    
        sum=sum+1;       
    end

    % Solve lambda    

    dfx=zeros(Nx,Ny);
    dfy=zeros(Nx,Ny);
    for i=1:Nx-1
        dfx(i,:)=(v(i+1,:)-v(i,:))/dx;
    end
    dfx(Nx,:)=dfx(Nx-1,:);
    for j=1:Ny-1
        dfy(:,j)=(v(:,j+1)-v(:,j))/dy;
    end
    dfy(:,Ny)=dfy(:,Ny-1);

    
   num=0;
   den=0;
   for i=1:Nx
        for j=1:Ny
            num=num+(dfx(i,j)^2+dfy(i,j)^2)*dx*dy;
            den=den+v(i,j)^3*dx*dy;
        end
   end
   lambda=num/den;


   ddf=zeros(Nx,Ny);
    for i=2:Nx-1
        for j=2:Ny-1
            ddf(i,j)=(v(i,j+1)+v(i,j-1)+v(i-1,j)+v(i+1,j)-4*v(i,j))/dy/dx;
        end
    end
   error=zeros(Nx,Ny);
    for i=1:Nx
        for j=1:Ny
            r_g=sqrt(x(i)^2+y(j)^2);
            if r_g>1
                continue
            end
            lhs=ddf(i,j);
            rhs=v(i,j)^2;
            error(i,j)=(lhs-rhs)/v(i,j);
        end
    end

   f=lambda*v;
   figure(1)
   surf(x,y,v)


end

aa=0;
E=0;
for i=1:Nx
    for j=1:Nx
        if error(i,j)<10
            aa=aa+1;
            E=E+error(i,j);
        end
    end
end
E=E/aa

