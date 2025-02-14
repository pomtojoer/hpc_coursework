clc;
close all;

%% Setting the graph interpreter as latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter','latex');

%% Defining the file names and locations to be plotted
re100 = "Data/1_1_161_161_100_data.txt";
re400 = "Data/1_1_161_161_400_data.txt";
re1000 = "Data/1_1_161_161_1000_data.txt";
re3200 = "Data/1_1_161_161_3200_data.txt";
lx1ly2 = "Data/1_2_161_161_100_data.txt";
lx2ly1 = "Data/2_1_161_161_100_data.txt";

filenames = [re100, "Re=100";
             re400, "Re=400";
             re1000, "Re=1000";
             re3200, "Re=3200"];

%% Plotting the required plots
xmin = zeros(4,1);
ymin = zeros(4,1);
         
for i=1:length(filenames)
    % opening the required files
    filename = filenames(i,1);
    plotlabel = filenames(i,2);
    fid = fopen(filename);
    
    % reading in the headers and parameters
    headers = fscanf(fid,"%s \n ",[1,7]);
    headersData = fscanf(fid,"%f, ",[1,7]);
    Lx = headersData(1);
    Ly = headersData(2);
    Nx = headersData(3);
    Ny = headersData(4);
    Re = headersData(7);
    
    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    
    % reading in the vorticity and streamfunctions matrices
    fscanf(fid,"%*s",1);
    w = fscanf(fid,"%f",[Ny,Nx])';  % vorticity
    fscanf(fid,"%*s",1);
    s = fscanf(fid,"%f",[Ny,Nx])';  % streamfunction
    
    % calculating the vertical velocity v = (s(i+1,j)-s(i-1,j))/(2*dx)
    if mod(Ny,2) == 1
        % odd number of grid points
        midpt_y = (Ny-1)/2+1;
    else
        % even number of grid points
        midpt_y = Ny/2;
    end
    v = -(s(midpt_y,3:Nx)-s(midpt_y,1:Nx-2))/2/dx;
    
    % calculating the vertical velocity u = (s(i,j+1)-s(i,j-1))/(2*dy)
    if mod(Nx,2) == 1
        % odd number of grid points
        midpt_x = (Nx-1)/2+1;
    else
        % even number of grid points
        midpt_x = Nx/2;
    end
    u = (s(3:Ny,midpt_x)-s(1:Ny-2,midpt_x))/2/dy;

    % determining the grid spacing
    xgrid = linspace(0,Lx,Nx);
    ygrid = linspace(0,Ly,Ny);
    
    % plotting u vs y
    figure(1)
    hold on
    plot(ygrid(2:Ny-1),u);
    
    % plotting v vs x
    figure(2)
    hold on
    plot(xgrid(2:Ny-1),v);
    
    % plotting vorticity and streamfunction contour for Re=100
    if i==1
        % generating mesh of xy location in grid
        [xcoord,ycoord] = meshgrid(xgrid,ygrid);
        
        figure(3)   % vorticity
        contourf(xcoord,ycoord,w,100);
        
        figure(4)   % streamfunction
        contourf(xcoord,ycoord,s);
    end
    
    % finding the xy coord of the minimum streamfunction
    xmin(i) = xcoord(s==min(s,[],'all'));
    ymin(i) = ycoord(s==min(s,[],'all'));
    
    % plotting the vorticity and streamfunction contours in the same plot
    % for comparison with one another
    splot = str2num("24" + (i*2));
    wplot = str2num("24" + (i*2-1));
    figure(10)
    hold on
    subplot(splot), contourf(s), axis('square'), title("streamfunction @"+plotlabel);
    subplot(wplot), contourf(w,100), axis('square'), title("vorticity @"+plotlabel);    
end

%% Adding title, legend and axis to plots then saving them
figure(1)
title('Plot of $u$ agaisnt $y$ for $x=0.5,\,L_x=1,\,L_y=1$ on a 161x161 grid');
grid on
grid minor
xlabel("$y$");
ylabel("$u$");
legend(filenames(:,2), "Location", "southwest");
% saveas(gcf,'Images/u_vs_y.png');

figure(2)
title('Plot of $v$ agaisnt $x$ for $y=0.5,\,L_x=1,\,L_y=1$ on a 161x161 grid');
grid on
grid minor
xlabel("$x$");
ylabel("$v$");
legend(filenames(:,2), "Location", "northwest");
% saveas(gcf,'Images/v_vs_x.png');

figure(3)
title('Contour plot for vorticity at $Re=100,\,L_x=1,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
% colorbar("eastoutside");
saveas(gcf,'Images/w_100.png');

figure(4)
title('Contour plot for streamfunction at $Re=100,\,L_x=1,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
% saveas(gcf,'Images/s_100.png');

figure(10)
% saveas(gcf,'Images/collated.png');

%% Plotting the plots with different Lx and Ly
filenames2 = [lx1ly2, lx2ly1];

for i=1:length(filenames2)
    % opening the required files
    filename = filenames2(i);
    fid = fopen(filename);
    
    % reading in the headers and parameters
    headers = fscanf(fid,"%s \n ",[1,7]);
    headersData = fscanf(fid,"%f, ",[1,7]);
    Lx = headersData(1);
    Ly = headersData(2);
    Nx = headersData(3);
    Ny = headersData(4);
    Re = headersData(7);
    
    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    
    % reading in the vorticity and streamfunctions matrices
    fscanf(fid,"%*s",1);
    w = fscanf(fid,"%f",[Ny,Nx])';  % vorticity
    fscanf(fid,"%*s",1);
    s = fscanf(fid,"%f",[Ny,Nx])';  % streamfunction
    
    % generating mesh of xy location in grid
    xgrid = linspace(0,Lx,Nx);
    ygrid = linspace(0,Ly,Ny);
    [xcoord,ycoord] = meshgrid(xgrid,ygrid);
    
    % plotting the streamfunction plots
    figure(5+i)
    contourf(xcoord,ycoord,s);
    daspect([1 1 1])
end

%% Adding title, legend and axis to plots then saving them
figure(6)
title('Contour plot for vorticity at $Re=100,\,L_x=1,\,L_y=2$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
% saveas(gcf,'Images/lx1ly2.png');

figure(7)
title('Contour plot for streamfunction at $Re=100,\,L_x=2,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
% saveas(gcf,'Images/lx2ly1.png');

%% Displaying the calculated streamfunction min coordinates
disp('x location of min');
disp(xmin');
disp('y location of min');
disp(ymin');

