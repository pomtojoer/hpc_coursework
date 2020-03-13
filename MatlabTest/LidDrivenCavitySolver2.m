clear ALL;
clc;
l=1; % Sides length of square
nx=5; %Number of grid in x-direction
ny=nx; % Number of grid in y-direction
nu=1/100; % 1/Reynolds number
dx=l/nx; % Element Size in x direction
dy=l/ny; % Element Size in y direction
U=1; % Top Velocity
rho=1.22; % Density
% Defining variables
si=zeros(nx,ny);
w=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
p=zeros(nx,ny);
for it=1 : 1000
%%%%%% Stream Function calculation using point gauss seidal %%%%%%%%%%%%%%%
for k=1:1000  
    dummy=si; % dummy variable used for convergence criteria
    for j=2:nx-1
        for i=2:ny-1
            si(j,i)=(0.25)*(((dx^2)*w(j,i))+ si(j,i+1)+ si(j,i-1)+si(j+1,i)+si(j-1,i));
        end
    end
    if abs(dummy-si)<=0.001 %%%Convergence criteria
    break;
    end
end
%%%%%%%%%%%% vorticity calculation using point gauss seidal method %%%%%%%%
dummy2=w;
for j=2:nx-1
    for i=2:ny-1
        w(j,i) = ((0.25)*(w(j,i+1)+w(j,i-1)+w(j+1,i)+w(j-1,i)))  - (  (0.25*dx/nu)*((si(j+1,i)-si(j-1,i))/2)*((w(j,i+1)-w(j,i-1))/2)  )   +  (  (0.25*dx/nu)*((si(j,i+1)-si(j,i-1))/2)*((w(j+1,i)-w(j-1,i))/2)  );
    end
end
%%%%% Boundary condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w(ny,:)=((2/(dx*dx))*(si(ny,:)-si(ny-1,:)))-(2*U/dx);
w(1,:)=(2/(dx*dx))*(si(1,:)-si(2,:));
w(:,nx)=(2/(dx*dx))*(si(:,nx)-si(:,nx-1));
w(:,1)=(2/(dx*dx))*(si(:,1)-si(:,2));
si(:,nx)=0;
si(:,1)=0;
si(1,:)=0;
si(ny,:)=0;
if abs(dummy2-w)<=0.001 %%% vortex convergence criteria
    break;
end
end
u(1,:)=0;
v(1,:)=0;
u(ny,:)=1;
v(ny,:)=0;
u(:,1)=0;
v(:,1)=0;
u(:,nx)=0;
v(:,nx)=0;
for j=2:nx-1
    for i=2:ny-1
    u(j,i)=(si(j+1,i)-si(j-1,i))/(2*dy);
    v(j,i)=-(si(j,i+1)-si(j,i-1))/(2*dx);
    end
end
for k=1:1000
dummyp=p;
for j=2:nx-1
    for i=2:ny-1
p(j,i)=((0.25)*(p(j,i+1)+p(j,i-1)+ p(j+1,i)+p(j-1,i))) - ((rho/(2*dx*dx))*(si(j,i+1)-2*si(j,i)+si(j,i-1))*(si(j+1,i)-2*si(j,i)+si(j-1,i))) + ( ((rho)/(32*dx*dx))*((si(j+1,i+1)-si(j-1,i+1)-si(j+1,i-1)+si(j-1,i-1))^2)) ;
    end
end
if abs(dummyp-p)<=0.001 %%%Convergence criteria
    break;
end
end
subplot(2,3,1);
plot(linspace(0,1,nx),v((ny/2)+0.5,:));
xlabel('X-direction');
ylabel('v-velocity');
title('X-direction vs V velocity')
subplot(2,3,2);
plot(u(:,(nx/2)+0.5),(linspace(0,1,ny))');
xlabel('U-velocity');
ylabel('Y-direction');
title(' U velocity vs Y-direction')
subplot(2,3,3), contourf(w,20);
xlabel('nx'),ylabel('ny'),title('Vorticity');axis('square','tight');colorbar
title('Vorticity')
subplot(2,3,4), contourf(si);
xlabel('nx'),ylabel('ny'),title('Stream function');axis('square','tight');colorbar
title('stream function')
subplot(2,3,5), contourf(p);
xlabel('nx'),ylabel('ny'),title('Pressure');axis('square','tight');colorbar
title('pressure') 
