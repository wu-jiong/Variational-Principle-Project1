clear;clc

h=0.01;     
a=1;      
b=1;     
m=a/h;    
n=b/h;     
f=zeros(n+1,m+1); 


x = 0:h:a;
y = 0:h:b;
Nx=length(x);
Ny=length(y);
dx=h;
dy=h;

for i=1:Nx
    for j=1:Ny
        f(i,j)=(x(i)-0.5).^2+(y(j)-0.5).^2;
    end
end
sigma=0.1;
f=-f/(2*sigma^2);
f=exp(f)/(2*pi*sigma^2);
surf(x,y,f)
title('Iterations k=0',FontSize=14)
lambda=1;

q=1.2;
for kk=1:10
    s=1; 
    sum=0; 
    v=zeros(n+1,m+1);
    while s>1e-6   
        w=v;    
        for k=2:1:n
            for l=2:1:m
                v(k,l)=(1-q)*v(k,l)+q*(v(k-1,l)+v(k,l-1)+v(k+1,l)+v(k,l+1)-h*h*(f(k,l)-tanh(f(k,l))^2))/4;
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
            num=num+(v(i,j)^2+dfx(i,j)^2+dfy(i,j)^2)*dx*dy;
            den=den+v(i,j)^3*dx*dy;
        end
   end
   lambda=num/den;


    real_lambda=calculate_lambda(v,dfx,dfy,x,y,2);
    error=check_error(v,x,y,2);

    f=lambda*v;
    surf(x,y,-v)
    title(['Interations k=',num2str(kk)],FontSize=14)

    a=0;
    E=0;
    for i=2:Nx-1
        for j=2:Nx-1
                a=a+1;
                E=E+abs(error(i,j));
        end
    end
    E=E/a;

end




