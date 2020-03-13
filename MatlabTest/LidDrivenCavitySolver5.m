clc;
clear; 

Lx = 1;
Ly = 1;
Nx = 5;
Ny = 5;
dt = 0.0001;
T = 100.0;
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

psiVarInterior = zeros(interiorNy,interiorNx);
psiVarLeft = zeros(1,interiorNy)';
psiVarRight = zeros(1,interiorNy)';
psiVarTop = zeros(1,interiorNx);
psiVarBottom = zeros(1,interiorNx);

psiVarInterior = [7     78    87  ;  
49    23    3    ;
73    9    27 ];
psiVarLeft = [30,65,40];
psiVarRight = [58,40,29];
psiVarBottom = [44,42,3];
psiVarTop = [72,92,12];

psiVarLeft = psiVarLeft';
psiVarRight = psiVarRight';

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

totalT = 0;
while totalT < T
    % Applying BC
    omegaTop = (psiVarTop-psiVarInterior(interiorNy,:))*2/dy/dy-2*U/dy;   % Top
    omegaBottom = (psiVarBottom-psiVarInterior(1,:))*2/dy/dy;       % Bottom
    omegaLeft = (psiVarLeft-psiVarInterior(:,1))*2/dx/dx;       % Left
    omegaRight = (psiVarRight-psiVarInterior(:,interiorNy))*2/dx/dx;      % Right
    
    % calculating interior vorticity at t
    interiorpsiVarV = psiVarInterior(:);
    temp = LapacianMatrix*interiorpsiVarV;
    temp = reshape(temp,[interiorNx,interiorNy]);
    temp(:,1) = temp(:,1) - psiVarLeft/dx/dx; % left interior
    temp(1,:) = temp(1,:) - psiVarBottom/dy/dy; % bottom interior
    temp(:,interiorNx) = temp(:,interiorNx) - psiVarRight/dx/dx; % right interior
    temp(interiorNy,:) = temp(interiorNy,:) - psiVarTop/dy/dy; % top interior
    omegaInterior = temp;
    
    % calculating interior vorticity at t+dt
    tempMat1 = zeros(interiorNy,interiorNx);
    tempMat2 = zeros(interiorNy,interiorNx);
    tempMat3 = zeros(interiorNy,interiorNx);
    tempMat4 = zeros(interiorNy,interiorNx);
    
    tempMat1(:,2:interiorNx-1) = psiVarInterior(:,3:interiorNx)-psiVarInterior(:,1:interiorNx-2);
    tempMat1(:,1) = psiVarInterior(:,2)-psiVarLeft;
    tempMat1(:,interiorNx) = psiVarRight-psiVarInterior(:,interiorNx-1);
    
    tempMat2(2:interiorNy-1,:) = omegaInterior(3:interiorNy,:)-omegaInterior(1:interiorNy-2,:);
    tempMat2(1,:) = omegaInterior(2,:)-omegaBottom;
    tempMat2(interiorNy,:) = omegaTop-omegaInterior(interiorNy-1,:);
    
    tempMat3(2:interiorNy-1,:) = psiVarInterior(3:interiorNy,:)-psiVarInterior(1:interiorNy-2,:);
    tempMat3(1,:) = psiVarInterior(2,:)-psiVarBottom;
    tempMat3(interiorNy,:) = psiVarTop-psiVarInterior(interiorNy-1,:);
    
    tempMat4(:,2:interiorNx-1) = omegaInterior(:,3:interiorNx)-omegaInterior(:,1:interiorNx-2);
    tempMat4(:,1) = omegaInterior(:,2)-omegaLeft;
    tempMat4(:,interiorNx) = omegaRight-omegaInterior(:,interiorNx-1);

    mat = (tempMat1.*tempMat2 - tempMat3.*tempMat4)*dt/4/dx/dy;
    interiorOmegaV = omegaInterior(:);
    temp = LapacianMatrix*interiorOmegaV;
    temp = reshape(temp,[interiorNx,interiorNy]);
    temp(1:interiorNy,1) = temp(1:interiorNy,1) - omegaLeft/dx/dx; % left interior
    temp(1:interiorNy,interiorNx) = temp(1:interiorNy,interiorNx) - omegaRight/dx/dx; % right interior
    temp(1,1:interiorNx) = temp(1,1:interiorNx) - omegaBottom/dy/dy; % bottom interior
    temp(interiorNy,1:interiorNx) = temp(interiorNy,1:interiorNx) - omegaTop/dy/dy; % top interior
    omegaInterior = -temp*dt/Re + mat + omegaInterior;
    
    % setting up poisson
    temp = zeros(interiorNy,interiorNx);
    temp(:,1) = omegaInterior(:,1) + psiVarLeft/dx/dx; % left interior
    temp(1,:) = omegaInterior(1,:) + psiVarBottom/dy/dy; % bottom interior
    temp(:,interiorNx) = omegaInterior(:,interiorNx) + psiVarRight/dx/dx; % right interior
    temp(interiorNy,:) = omegaInterior(interiorNy,:) + psiVarTop/dy/dy; % top interior
    interiorOmegaV = temp(:);
    
    % solving poisson
    t = linsolve(LapacianMatrix,interiorOmegaV);
    psiVarInterior = reshape(t,[interiorNx,interiorNy]);
    
    combinedOmega = zeros(Ny,Nx);
    combinedOmega(Ny,2:Nx-1) = omegaTop;
    combinedOmega(1,2:Nx-1) = omegaBottom;
    combinedOmega(2:Ny-1,1) = omegaLeft;
    combinedOmega(2:Ny-1,Nx) = omegaRight;
    combinedOmega(2:Ny-1,2:Nx-1) = omegaInterior;
    
    combinedPsiVar = zeros(Ny,Nx);
    combinedPsiVar(Ny,2:Nx-1) = psiVarTop;
    combinedPsiVar(1,2:Nx-1) = psiVarBottom;
    combinedPsiVar(2:Ny-1,1) = psiVarLeft;
    combinedPsiVar(2:Ny-1,Nx) = psiVarRight;
    combinedPsiVar(2:Ny-1,2:Nx-1) = psiVarInterior;
    
    psiVarTop
    psiVarBottom
    psiVarLeft
    psiVarRight
    psiVarInterior
    omegaTop
    omegaBottom
    omegaLeft
    omegaRight
    omegaInterior
    
    subplot(121), contour((fliplr(combinedOmega))), axis('square');            % Vorticity
    subplot(122), contour((fliplr(combinedPsiVar))), axis('square');          % Stream function
    pause(0.01)
    
    totalT = totalT + dt;
end




% omegaMat(2:Ny-1,2:Nx-1) = temp;

% solving stream function
