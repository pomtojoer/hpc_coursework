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

psiInterior = zeros(interiorNy,interiorNx);
psiLeft = zeros(1,interiorNy)';
psiRight = zeros(1,interiorNy)';
psiTop = zeros(1,interiorNx);
psiBottom = zeros(1,interiorNx);
% 
% psiInterior = [7     78    87  ;  
% 49    23    3    ;
% 73    9    27 ];
% psiLeft = [30,65,40];
% psiRight = [58,40,29];
% psiBottom = [44,42,3];
% psiTop = [72,92,12];
% 
% psiLeft = psiLeft'
% psiRight = psiRight'

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

LapacianMatrix2 = zeros((Ny)*(Ny),(Nx)*(Nx));
for j = 1:(Nx)*(Nx)
    for i = 1:(Ny)*(Ny)
        if i==j
            LapacianMatrix2(i,j) = gamma;
        elseif (i==j-(Ny)) || (i-(Ny)==j)
            LapacianMatrix2(i,j) = -alpha;
        elseif (((i==j-1) && (mod(i,Ny)~=0)) || ((i-1==j) && (mod(j,Ny)~=0)))
            LapacianMatrix2(i,j) = -beta;
        end
    end
end

l=1; % Sides length of square
nx=5; %Number of grid in x-direction
ny=nx; % Number of grid in y-direction
nu=1/100; % 1/Reynolds number
dx=l/(nx-1); % Element Size in x direction
dy=l/(ny-1); % Element Size in y direction
U=1; % Top Velocity
rho=1.22; % Density
% Defining variables
si=zeros(nx,ny);
w=zeros(nx,ny);

totalT = 0;
while totalT < T
    disp(totalT);
    w
    si
    % Applying BC
    omegaTop = (psiTop-psiInterior(interiorNy,:))*2/dy/dy-2*U/dy;   % Top
    omegaBottom = (psiBottom-psiInterior(1,:))*2/dy/dy;       % Bottom
    omegaLeft = (psiLeft-psiInterior(:,1))*2/dx/dx;       % Left
    omegaRight = (psiRight-psiInterior(:,interiorNy))*2/dx/dx;      % Right
    
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
    temp(:,1) = temp(:,1) - psiLeft/dx/dx; % left interior
    temp(1,:) = temp(1,:) - psiBottom/dy/dy; % bottom interior
    temp(:,interiorNx) = temp(:,interiorNx) - psiRight/dx/dx; % right interior
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

    omega;
    
    me = zeros(Ny,Nx);
    me(Ny,2:Nx-1) = omegaTop;
    me(1,2:Nx-1) = omegaBottom;
    me(2:Ny-1,1) = omegaLeft;
    me(2:Ny-1,Nx) = omegaRight;
    me(2:Ny-1,2:Nx-1) = omegaInterior;
    
    omegaInterior
    % temp = LapacianMatrix*interiorPsiV;
    % temp = reshape(temp,[interiorNx,interiorNy]);
    omegaInterior(:,1) = omegaInterior(:,1) + psiLeft/dx/dx; % left interior
    omegaInterior(1,:) = omegaInterior(1,:) + psiBottom/dy/dy; % bottom interior
    omegaInterior(:,interiorNx) = omegaInterior(:,interiorNx) + psiRight/dx/dx; % right interior
    omegaInterior(interiorNy,:) = omegaInterior(interiorNy,:) + psiTop/dy/dy; % top interior
    interiorOmegaV = omegaInterior(:);
    omegaInterior
    
    t = linsolve(LapacianMatrix,interiorOmegaV);
    psiInterior = reshape(t,[interiorNx,interiorNy]);
    
    me2 = zeros(Ny,Nx);
    me2(Ny,2:Nx-1) = psiTop;
    me2(1,2:Nx-1) = psiBottom;
    me2(2:Ny-1,1) = psiLeft;
    me2(2:Ny-1,Nx) = psiRight;
    me2(2:Ny-1,2:Nx-1) = psiInterior;
    
    meV = me(:)
    t2 = linsolve(LapacianMatrix2,meV);
    psiInterior2 = reshape(t2,[Nx,Ny]);
    
    subplot(121), contour((fliplr(me))), axis('square');            % Vorticity
    subplot(122), contour((fliplr(me2))), axis('square');          % Stream function
    pause(0.01)
    psiInterior;
    
    totalT = totalT + dt;
    
    %%%%%% Stream Function calculation using point gauss seidal %%%%%%%%%%%%%%%
    for k=1:1000  
        dummy=si; % dummy variable used for convergence criteria
        for j=2:nx-1
            for i=2:ny-1
                si(j,i)=(0.25)*(((dx^2)*w(j,i))+ si(j,i+1)+ si(j,i-1)+si(j+1,i)+si(j-1,i));
            end
        end
        if abs(dummy-si)<=0.001 %%%Convergence criteria
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
    
end




% omegaMat(2:Ny-1,2:Nx-1) = temp;

% solving stream function
