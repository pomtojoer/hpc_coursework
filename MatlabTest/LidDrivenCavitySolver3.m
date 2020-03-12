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

psi = [7     72    40    3     3   ;  
49    44    65    27    69    ;
73    78    92    29    9     ;
58    23    42    40    57    ;
30    9     87    12    60 ;];

alpha = 1/dx/dx;
beta = 1/dy/dy;
gamma = 2*(alpha+beta);

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
omegaMat = omega;

% calculating interior vorticity at t
for i = 2:Ny-1
    for j = 2:Nx-1
        omega(j,i) = -1*((psi(j,i+1)-2*psi(j,i)+psi(j,i-1))/dx/dx+(psi(j+1,i)-2*psi(j,i)+psi(j-1,i))/dy/dy);
    end
end

interiorPsi = psi(2:Ny-1,2:Nx-1);
interiorPsiV = interiorPsi(:);
temp = LapacianMatrix*interiorPsiV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp
temp(1:interiorNy,1) = temp(1:interiorNy,1) - psi(2:Ny-1,1)/dx/dx; % left interior
temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - psi(2:Ny-1,Nx)/dx/dx; % right interior
temp(1,1:interiorNx) = temp(1,1:interiorNx) - psi(1,2:Nx-1)/dy/dy; % bottom interior
temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - psi(Ny,2:Nx-1)/dy/dy; % right interior
omegaMat(2:Ny-1,2:Nx-1) = temp;

omega
omegaMat
% calculating interior vorticity at t+dt
term1Mat = zeros(Ny,Nx);
term2Mat = zeros(Ny,Nx);
omegaNew = zeros(Ny,Nx);

for j = 2:Nx-1
    for i = 2:Ny-1
        term1 = dt/4/dx/dy*((psi(j,i+1)-psi(j,i-1))*(omega(j+1,i)-omega(j-1,i))-(psi(j+1,i)-psi(j-1,i))*(omega(j,i+1)*omega(j,i-1)));
        term1Mat(j,i) = term1;
        term2 = (dt/Re)*((omega(j,i+1)-2*omega(j,i)+omega(j,i-1))/dx/dx + (omega(j+1,i)-2*omega(j,i)+omega(j-1,i))/dy/dy);
        term2Mat(j,i) = term2;
        omegaNew(i,j) = term1 + term2 + omega(i,j);
    end
end
term1Mat


tempMat1 = psi(2:Ny-1,3:Nx)-psi(2:Ny-1,1:Nx-2);
tempMat2 = omegaMat(3:Ny,2:Nx-1)-omegaMat(1:Ny-2,2:Nx-1);
tempMat3 = psi(3:Ny,2:Nx-1)-psi(1:Ny-2,2:Nx-1);
tempMat4 = omegaMat(2:Ny-1,3:Nx)-omegaMat(2:Ny-1,1:Nx-2);

mat = (tempMat1.*tempMat2 - tempMat3.*tempMat4)*dt/4/dx/dy

interiorOmega = omegaMat(2:Ny-1,2:Nx-1);
interiorOmegaV = interiorOmega(:);
temp = LapacianMatrix*interiorOmegaV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp(1:interiorNy,1) = temp(1:interiorNy,1) - omegaMat(2:Ny-1,1)/dx/dx; % left interior
temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - omegaMat(2:Ny-1,Nx)/dx/dx; % right interior
temp(1,1:interiorNx) = temp(1,1:interiorNx) - omegaMat(1,2:Nx-1)/dy/dy; % bottom interior
temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - omegaMat(Ny,2:Nx-1)/dy/dy; % right interior
temp = -temp*dt/Re;
% omegaMat(2:Ny-1,2:Nx-1) = temp;
temp
term2Mat
% solving stream function
