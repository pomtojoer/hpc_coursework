clc;
clear; 

%% Variables to be defined 
Lx = 1;         % domain length in the x direction
Ly = 1;         % domain length in the y direction
Nx = 11;        % grid size in the x direction
Ny = 5;         % grid size in the y direction
dt = 0.0005;    % time strp
T = 1.0;        % final time
Re = 100;       % Reynolds number

%% Checking the stability criterion
U = 1.0;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

if Re*dx*dy/4 <= dt
    disp("Choose smaller time step");
    return
end

%% Initialising the streamfunction and vorticity matrices
interiorNy = Ny-2;
interiorNx = Nx-2;

w = zeros(Ny,Nx);
s = zeros(Ny,Nx);

%% Generating the Laplacian matrix for solve
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

%% Solving the lid driven cavity problem
iterations = 0;
totalT = 0;
while totalT < T
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
    
    
    % ##################### animated plotting ##################### 
    figure(1)
    subplot(221), contour(((w))), axis('square');   % Vorticity
    subplot(222), contour(((s))), axis('square');	% Stream function
    pause(0.0001)
    
    totalT = totalT + dt;
    iterations = iterations + 1;
end

%% Displaying and plotting the final values
disp("The vorticity matrix looks like: ")
disp(w)
disp("The streamfunction matrix looks like: ")
disp(s)
figure(2)
subplot(221), contour(((w))), axis('square');	% Vorticity
subplot(222), contour(((s))), axis('square');	% Streamfunction