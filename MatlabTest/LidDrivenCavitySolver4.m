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


omegaInterior = zeros(interiorNy,interiorNx);
omegaLeft = zeros(interiorNy);
omegaRight = zeros(interiorNy);
omegaTop = zeros(interiorNx);
omegaBottom = zeros(interiorNx);

psiInterior = zeros(interiorNy,interiorNx);
psiLeft = zeros(1,interiorNy);
psiRight = zeros(1,interiorNy);
psiTop = zeros(1,interiorNx);
psiBottom = zeros(1,interiorNx);

psiInterior = [7     78    87  ;  
49    23    3    ;
73    9    27 ];
psiLeft = [30,65,40];
psiRight = [58,40,29];
psiBottom = [44,42,3];
psiTop = [72,92,12];

omega = zeros(Ny,Nx);
psi = zeros(Ny,Nx);
psi(Ny,2:Nx-1) = psiTop;
psi(1,2:Nx-1) = psiBottom;
psi(2:Ny-1,1) = psiLeft;
psi(2:Ny-1,Nx) = psiRight;
psi(2:Ny-1,2:Nx-1) = psiInterior;

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

% Applying BC
omegaTop = (psiTop-psiInterior(interiorNy,:))*2/dy/dy-2*U/dy;   % Top
omegaBottom = (psiBottom-psiInterior(1,:))*2/dy/dy;       % Bottom
omegaLeft = (psiLeft'-psiInterior(:,1))*2/dx/dx;       % Left
omegaRight = (psiRight'-psiInterior(:,interiorNy))*2/dx/dx;      % Right

omega(Ny,:) = (psi(Ny,:)-psi(Ny-1,:))*2/dy/dy-2*U/dy;   % Top
omega(1,:) = (psi(1,:)-psi(2,:))*2/dy/dy;       % Bottom
omega(:,1) = (psi(:,1)-psi(:,2))*2/dx/dx;       % Left
omega(:,Nx) = (psi(:,Nx)-psi(:,Nx-1))*2/dx/dx;      % Right

% calculating interior vorticity at t
for i = 2:Ny-1
    for j = 2:Nx-1
        omega(j,i) = -1*((psi(j,i+1)-2*psi(j,i)+psi(j,i-1))/dx/dx+(psi(j+1,i)-2*psi(j,i)+psi(j-1,i))/dy/dy);
    end
end

interiorPsiV = psiInterior(:);
temp = LapacianMatrix*interiorPsiV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp(:,1) = temp(:,1) - psiLeft'/dx/dx; % left interior
temp(1,:) = temp(1,:) - psiBottom/dy/dy; % bottom interior
temp(:,interiorNx) = temp(:,interiorNx) - psiRight'/dx/dx; % right interior
temp(interiorNy,:) = temp(interiorNy,:) - psiTop/dy/dy; % top interior
omegaInterior = temp;

% calculating interior vorticity at t+dt
tMat = zeros(Ny,Nx);
term1Mat = zeros(Ny,Nx);
term2Mat = zeros(Ny,Nx);
omegaNew = zeros(Ny,Nx);
omegaNew2 = zeros(Ny,Nx);

for j = 2:Nx-1
    for i = 2:Ny-1
        tMat(j,i) = (omega(j,i+1)-omega(j,i-1));
        term1 = dt/4/dx/dy*((psi(j,i+1)-psi(j,i-1))*(omega(j+1,i)-omega(j-1,i))-(psi(j+1,i)-psi(j-1,i))*(omega(j,i+1)-omega(j,i-1)));
        term1Mat(j,i) = term1;
        term2 = (dt/Re)*((omega(j,i+1)-2*omega(j,i)+omega(j,i-1))/dx/dx + (omega(j+1,i)-2*omega(j,i)+omega(j-1,i))/dy/dy);
        term2Mat(j,i) = term2;
        omegaNew(j,i) = term1 + term2 + omega(j,i);
    end
end

tempMat1 = psi(2:Ny-1,3:Nx)-psi(2:Ny-1,1:Nx-2);
tempMat2 = omega(3:Ny,2:Nx-1)-omega(1:Ny-2,2:Nx-1);
tempMat3 = psi(3:Ny,2:Nx-1)-psi(1:Ny-2,2:Nx-1);
tempMat4 = omega(2:Ny-1,3:Nx)-omega(2:Ny-1,1:Nx-2);

mat = (tempMat1.*tempMat2 - tempMat3.*tempMat4)*dt/4/dx/dy;
interiorOmegaV = omegaInterior(:);
temp = LapacianMatrix*interiorOmegaV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp(1:interiorNy,1) = temp(1:interiorNy,1) - omega(2:Ny-1,1)/dx/dx; % left interior
temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - omega(2:Ny-1,Nx)/dx/dx; % right interior
temp(1,1:interiorNx) = temp(1,1:interiorNx) - omega(1,2:Nx-1)/dy/dy; % bottom interior
temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - omega(Ny,2:Nx-1)/dy/dy; % right interior
omegaInterior = -temp*dt/Re + mat + omegaInterior;


interiorOmega = omega(2:Ny-1,2:Nx-1);
interiorOmegaV = interiorOmega(:);
temp = LapacianMatrix*interiorOmegaV;
temp = reshape(temp,[interiorNx,interiorNy]);
temp(1:interiorNy,1) = temp(1:interiorNy,1) - omega(2:Ny-1,1)/dx/dx; % left interior
temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - omega(2:Ny-1,Nx)/dx/dx; % right interior
temp(1,1:interiorNx) = temp(1,1:interiorNx) - omega(1,2:Nx-1)/dy/dy; % bottom interior
temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - omega(Ny,2:Nx-1)/dy/dy; % right interior
temp = -temp*dt/Re;
omega(2:Ny-1,2:Nx-1) = omega(2:Ny-1,2:Nx-1) + temp + mat;

omegaInterior

% omegaMat(2:Ny-1,2:Nx-1) = temp;

% solving stream function
