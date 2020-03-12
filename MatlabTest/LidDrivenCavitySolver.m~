clc;
clear; 

Lx = 1;
Ly = 1;
Nx = 5;
Ny = 5;
dt = 0.0001;
T = 1.0;
Re = 100;

U = 1.0;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);


interiorNx = Nx-2;
interiorNy = Ny-2;



omega = zeros(Ny,Nx);
psi = zeros(Ny,Nx);

counter = 1;
for i=1:Ny
    for j=1:Nx
        psi(i,j) = counter;
        counter = counter + 1;
    end
end

% Iterative method

% while T > 0
%     
%     
% end





% Applying BC
omega(Ny,:) = (psi(Ny,:)-psi(Ny-1,:))*2/dy/dy-2*U/dy;   % Top
omega(1,:) = (psi(1,:)-psi(2,:))*2/dy/dy;       % Bottom
omega(:,1) = (psi(:,1)-psi(:,2))*2/dx/dx;       % Left
omega(:,Nx) = (psi(:,Nx)-psi(:,Nx-1))*2/dx/dx;      % Right


% creating vector for matrix multiplication
temp = psi(2:Ny-1,2:Nx-1);
% temp = psi'
psiV = temp(:);
% psiV = psi(:);

alpha = 1/dx/dx;
beta = 1/dy/dy;
gamma = 2*(alpha+beta);

r2 = gamma*ones(Ny,1);
r = -beta*ones(Ny-1,1);
A = diag(r2,0)+diag(r(1:Ny-1),1)+diag(r(1:Ny-1),-1);
B = -alpha*eye(Ny,Nx);

LapacianMatrix = zeros((interiorNy)*(interiorNy),(interiorNx)*(interiorNx));
for j = 1:(interiorNx)*(interiorNx)
    for i = 1:(interiorNy)*(interiorNy)
        if i==j
            LapacianMatrix(i,j) = gamma;
        elseif (i==j-(interiorNy)) || (i-(interiorNx)==j)
            LapacianMatrix(i,j) = -alpha;
        elseif (((i==j-1) && (mod(i,interiorNy)~=0)) || ((i-1==j) && (mod(j,interiorNx)~=0)))
            LapacianMatrix(i,j) = -beta;
        end
    end
end

temp = LapacianMatrix*psiV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp(1:interiorNy,1) = temp(1:interiorNy,1) - psi(2:Ny-1,1)/dx/dx; % left interior
temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - psi(2:Ny-1,Nx)/dx/dx; % right interior
temp(1,1:interiorNx) = temp(1,1:interiorNx) - psi(1,2:Nx-1)/dy/dy; % bottom interior
temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - psi(Ny,2:Nx-1)/dy/dy; % right interior


% calculating interior vorticity at t
for i = 2:Ny-1
    for j = 2:Nx-1
        omega(i,j) = -1*((psi(i+1,j)-2*psi(i,j)+psi(i-1,j))/dx/dx+(psi(i,j+1)-2*psi(i,j)+psi(i,j-1))/dy/dy);
    end
end

omega

% calculating interior vorticity at t+dt
for j = 2:Nx-1
    for i = 2:Ny-1
        term1 = dt/4/dx/dy*((psi(i+1,j)-psi(i-1,j))*(omega(i,j+1)-omega(i,j-1))-(psi(i,j+1)-psi(i,j-1))*(omega(i+1,j)*omega(i-1,j)));
        term2 = dt/Re*((omega(i+1,j)-2*omega(i,j)+omega(i-1,j))/dx/dx + (omega(i,j+1)-2*omega(i,j)+omega(i,j-1))/dy/dy);
        omega(i,j) = term1 + term2 + omega(i,j);
    end
end

% solving stream function
