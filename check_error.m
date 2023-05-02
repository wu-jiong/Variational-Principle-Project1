function [error] = check_error(f0,x,y,m)
    Nx=length(x);
    Ny=length(y);
    dx=x(2)-x(1);
    dy=y(2)-y(1);
    ddf=zeros(Nx,Ny);
    for i=2:Nx-1
        for j=2:Ny-1
            ddf(i,j)=(f0(i,j+1)+f0(i,j-1)+f0(i-1,j)+f0(i+1,j)-4*f0(i,j))/dy/dx;
        end
    end
    
    error=zeros(Nx,Ny);
    if m==1
    for i=1:Nx
        for j=1:Ny
            r_g=sqrt(x(i)^2+y(j)^2);
            if r_g>1
                continue
            end
            error(i,j)=(ddf(i,j)-f0(i,j)+tanh(f0(i,j))^2)/f0(i,j);
        end
    end
    end

    if m==2
    for i=1:Nx
        for j=1:Ny
            error(i,j)=(ddf(i,j)-f0(i,j)+tanh(f0(i,j))^2)/f0(i,j);
        end
    end
    end

    if m==3
    for i=1:Nx
        for j=1:Ny
            r_g=sqrt(x(i)^2+y(j)^2);
            if r_g>1
                continue
            end
            error(i,j)=(ddf(i,j)-f0(i,j)^2+tanh(f0(i,j)));
        end
    end
    end


end

