% clc;
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

if Re*dx*dy/4 <= dt
    disp("Choose smaller time step");
    return
end

interiorNy = Ny-2;
interiorNx = Nx-2;

w = zeros(Ny,Nx);
s = zeros(Ny,Nx);

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

s = [7     72    40    3     3   ;  
49    44    65    27    69    ;
73    78    92    29    9     ;
58    23    42    40    57    ;
30    9     87    12    60 ;];

psiVarInterior = s(2:Ny-1,2:Nx-1);
psiVarLeft = s(2:Ny-1,1);
psiVarRight = s(2:Ny-1,Nx);
psiVarTop = s(Ny,2:Nx-1);
psiVarBottom = s(1,2:Nx-1);

alpha = 1/dx/dx;
beta = 1/dy/dy;
gamma = 2*(alpha+beta);

LapacianMatrix = zeros((interiorNx)*(interiorNy),(interiorNx)*(interiorNy));
for j = 1:(interiorNy)*(interiorNx)
    for i = 1:(interiorNy)*(interiorNx)
        if i==j
            LapacianMatrix(j,i) = gamma;
        elseif (j==i+interiorNy) || (i==j+interiorNy)
            LapacianMatrix(i,j) = -alpha;
        elseif (((i==j-1) && (mod(i,interiorNy)~=0)) || ((i-1==j) && (mod(j,interiorNy)~=0)))
            LapacianMatrix(i,j) = -beta;
        end
    end
end


totalT = 0;
while totalT < T
    disp("############################################################")
    disp("initia;");
    disp(w);
    disp(s);
    
    % ##################### Applying BC ##################### 
    % iteratively
    w(Ny,:) = (s(Ny,:) - s(Ny-1,:))*2/dy/dy-2*U/dy; % Top
    w(1,:) = (s(1,:) - s(2,:))*2/dy/dy; % Bottom
    w(:,1) = (s(:,1) - s(:,2))*2/dx/dx; % Left
    w(:,Nx) = (s(:,Nx) - s(:,Nx-1))*2/dx/dx; % Right
    w(1,1) = 0;
    w(Ny,1) = 0;
    w(Ny,Nx) = 0;
    w(1,Ny) = 0;
    
    % matrix
    omegaTop = (psiVarTop-psiVarInterior(interiorNy,:))*2/dy/dy-2*U/dy;   % Top
    omegaBottom = (psiVarBottom-psiVarInterior(1,:))*2/dy/dy;       % Bottom
    omegaLeft = (psiVarLeft-psiVarInterior(:,1))*2/dx/dx;       % Left
    omegaRight = (psiVarRight-psiVarInterior(:,interiorNy))*2/dx/dx;      % Right
    
    % difference
    diffTop = norm(omegaTop - w(Ny,2:Nx-1));
    diffBottom = norm(omegaBottom - w(1,2:Nx-1));
    diffLeft = norm(omegaLeft - w(2:Ny-1,1));
    diffRight = norm(omegaRight - w(2:Ny-1,Nx));
    diffInterior = norm(omegaInterior - w(2:Ny-1,2:Nx-1));
    
    disp("after set bc");
    disp(w);
    disp(s);
    
    % ##################### calculating interior vorticity at t ##################### 
    % iteratively
    for i = 2:Nx-1
        for j = 2:Ny-1
            w(j,i) = -(s(j,i+1)-2*s(j,i)+s(j,i-1))/dx/dx -(s(j+1,i)-2*s(j,i)+s(j-1,i))/dy/dy;
        end
    end
    
    % matrix
    interiorpsiVarV = psiVarInterior(:);
    temp = LapacianMatrix*interiorpsiVarV;
    temp = reshape(temp,[interiorNx,interiorNy]);
    temp(:,1) = temp(:,1) - psiVarLeft/dx/dx; % left interior
    temp(1,:) = temp(1,:) - psiVarBottom/dy/dy; % bottom interior
    temp(:,interiorNx) = temp(:,interiorNx) - psiVarRight/dx/dx; % right interior
    temp(interiorNy,:) = temp(interiorNy,:) - psiVarTop/dy/dy; % top interior
    omegaInterior = temp;
    
    % difference
    diffTop = norm(omegaTop - w(Ny,2:Nx-1));
    diffBottom = norm(omegaBottom - w(1,2:Nx-1));
    diffLeft = norm(omegaLeft - w(2:Ny-1,1));
    diffRight = norm(omegaRight - w(2:Ny-1,Nx));
    diffInterior = norm(omegaInterior - w(2:Ny-1,2:Nx-1));
    
    disp("after set interior");
    disp(w);
    disp(s);
    
    % ##################### calculating interior vorticity at t+dt ##################### 
    % iteratively
    temp = zeros(Ny,Nx);
    tempMat1 = zeros(Ny,Nx);
    for i = 2:Nx-1
        for j = 2:Ny-1
            temp(j,i) = dt/4/dx/dy*((s(j,i+1)-s(j,i-1))*(w(j+1,i)-w(j-1,i))-...
                (s(j+1,i)-s(j-1,i))*(w(j,i+1)-w(j,i-1))) + ...
                dt/Re*((w(j,i+1)-2*w(j,i)+w(j,i-1))/dx/dx+(w(j+1,i)-2*w(j,i)+w(j-1,i))/dy/dy) + ...
                w(j,i);
            tempMat1(j,i) = (s(j+1,i)-s(j-1,i))*(w(j,i+1)-w(j,i-1));
        end
    end
    w(2:Ny-1,2:Nx-1) = temp(2:Ny-1,2:Nx-1);
    
    % matrix
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
    omegaInterior = -1*temp*dt/Re + mat + omegaInterior;
    
    % difference
    diffTop = norm(omegaTop - w(Ny,2:Nx-1));
    diffBottom = norm(omegaBottom - w(1,2:Nx-1));
    diffLeft = norm(omegaLeft - w(2:Ny-1,1));
    diffRight = norm(omegaRight - w(2:Ny-1,Nx));
    diffInterior = norm(omegaInterior - w(2:Ny-1,2:Nx-1));
    
    disp("after update");
    disp(w);
    disp(s);
    
    % ##################### solving poisson ##################### 
    % iteratively
    tempw = w(2:Ny-1,2:Nx-1);
    tempw(:,1) = tempw(:,1) + s(2:Ny-1,1)/dx/dx; % left
    tempw(1,:) = tempw(1,:) + s(1,2:Nx-1)/dy/dy; % bottom
    tempw(:,interiorNx) = tempw(:,interiorNx) + s(2:Ny-1,Nx)/dx/dx; % right
    tempw(interiorNy,:) = tempw(interiorNy,:) + s(Ny,2:Nx-1)/dy/dy; % top
    interiorOmegaV = tempw;
    t = linsolve(LapacianMatrix,interiorOmegaV(:));
    s(2:Ny-1,2:Nx-1) = reshape(t,[interiorNx,interiorNy]);
    
    % matrix
    temp = omegaInterior;
    temp(:,1) = temp(:,1) + psiVarLeft/dx/dx; % left interior
    temp(1,:) = temp(1,:) + psiVarBottom/dy/dy; % bottom interior
    temp(:,interiorNx) = temp(:,interiorNx) + psiVarRight/dx/dx; % right interior
    temp(interiorNy,:) = temp(interiorNy,:) + psiVarTop/dy/dy; % top interior
    interiorOmegaV = temp(:);
    
    t = linsolve(LapacianMatrix,interiorOmegaV);
    psiVarInterior = reshape(t,[interiorNx,interiorNy]);
    
    omega = zeros(Ny,Nx);
    omega(Ny,2:Nx-1) = omegaTop; % Top
    omega(1,2:Nx-1) = omegaBottom; % Bottom
    omega(2:Ny-1,1) = omegaLeft; % Left
    omega(2:Ny-1,Nx) = omegaRight; % Right
    omega(2:Ny-1,2:Nx-1) = omegaInterior; % Interior
    
    psiVar = zeros(Ny,Nx);
    psiVar(Ny,2:Nx-1) = psiVarTop; % Top
    psiVar(1,2:Nx-1) = psiVarBottom; % Bottom
    psiVar(2:Ny-1,1) = psiVarLeft; % Left
    psiVar(2:Ny-1,Nx) = psiVarRight; % Right
    psiVar(2:Ny-1,2:Nx-1) = psiVarInterior; % Interior
    
    % difference
    diffTop = norm(psiVarTop - s(Ny,2:Nx-1));
    diffBottom = norm(psiVarBottom - s(1,2:Nx-1));
    diffLeft = norm(psiVarLeft - s(2:Ny-1,1));
    diffRight = norm(psiVarRight - s(2:Ny-1,Nx));
    diffInterior = norm(psiVarInterior - s(2:Ny-1,2:Nx-1));
    
    
    % ##################### plotting ##################### 
    subplot(221), contour(((w))), axis('square');            % Vorticity
    subplot(222), contour(((s))), axis('square');          % Stream function

    subplot(223), contour(((omega))), axis('square');            % Vorticity
    subplot(224), contour(((psiVar))), axis('square');          % Stream function
    pause(0.0001)
    
    totalT = totalT + dt;
    disp("after solve");
    disp(w);
    disp(s);
end


% omegaMat(2:Ny-1,2:Nx-1) = temp;

% solving stream function
