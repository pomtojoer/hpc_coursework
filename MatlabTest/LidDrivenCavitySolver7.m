clc;
clear; 

Lx = 1;
Ly = 1;
Nx = 11;
Ny = 5;
dt = 0.0001;
% dt = 0.01;
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

% s = [7     72    40    3     3   ;  
% 49    44    65    27    69    ;
% 73    78    92    29    9     ;
% 58    23    42    40    57    ;
% 30    9     87    12    60 ;];

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

iterations = 0;
totalT = 0;
while totalT < T
%     disp(w)
    
    % ##################### calculating interior vorticity at t ##################### 
    % iteratively
    for i = 2:Nx-1
        for j = 2:Ny-1
            w(j,i) = -(s(j,i+1)-2*s(j,i)+s(j,i-1))/dx/dx -(s(j+1,i)-2*s(j,i)+s(j-1,i))/dy/dy;
        end
    end
    
    % ##################### Applying BC ##################### 
    % iteratively
    w(Ny,:) = (s(Ny,:) - s(Ny-1,:))*2/dy/dy-2*U/dy; % Top
    w(1,:) = (s(1,:) - s(2,:))*2/dy/dy; % Bottom
    w(:,1) = (s(:,1) - s(:,2))*2/dx/dx; % Left
    w(:,Nx) = (s(:,Nx) - s(:,Nx-1))*2/dx/dx; % Right
    w(1,1) = 0;
    w(Ny,1) = 0;
    w(Ny,Nx) = 0;
    w(1,Nx) = 0;
    
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
    
    % ##################### solving poisson ##################### 
    % iteratively
    interiorOmegaV = w(2:Ny-1,2:Nx-1);
    t = linsolve(LapacianMatrix,interiorOmegaV(:));
    s(2:Ny-1,2:Nx-1) = reshape(t,[interiorNy,interiorNx]);
    
    
%     % ##################### plotting ##################### 
%     subplot(221), contour(((w))), axis('square');            % Vorticity
%     subplot(222), contour(((s))), axis('square');          % Stream function
%     pause(0.0001)
%     
    totalT = totalT + dt;
    iterations = iterations + 1;
end

w
s
subplot(221), contour(((w))), axis('square');	% Vorticity
subplot(222), contour(((s))), axis('square');	% Streamfunction

% omegaMat(2:Ny-1,2:Nx-1) = temp;

% solving stream function
