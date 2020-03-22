clc;
close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter','latex');

test = "1_1_5_5_100_data.txt";
re100 = "1_1_161_161_100_data.txt";
re400 = "1_1_161_161_400_data.txt";
re1000 = "1_1_161_161_1000_data.txt";
re3200 = "1_1_161_161_3200_data.txt";
lx1ly2 = "1_2_161_161_100_data.txt";
lx2ly1 = "2_1_161_161_100_data.txt";

filenames = [re100, "Re=100";
             re400, "Re=400";
             re1000, "Re=1000";
             re3200, "Re=3200"];

xmin = zeros(4,1);
ymin = zeros(4,1);
         
for i=1:length(filenames)
    filename = filenames(i,1);
    plotlabel = filenames(i,2);
    fid = fopen(filename);
    headers = fscanf(fid,"%s \n ",[1,7]);
    headersData = fscanf(fid,"%f, ",[1,7]);
    Lx = headersData(1);
    Ly = headersData(2);
    Nx = headersData(3);
    Ny = headersData(4);
    Re = headersData(7);
    
    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    
    fscanf(fid,"%*s",1);
    w = fscanf(fid,"%f",[Ny,Nx])';
    fscanf(fid,"%*s",1);
    s = fscanf(fid,"%f",[Ny,Nx])';
    
    if mod(Nx,2) == 1
        % odd number of grid points
        lhs_idx = (Nx-1)/2;
        rhs_idx = lhs_idx + 2;
    else
        % even number of grid points
        lhs_idx = Nx/2;
        rhs_idx = lhs_idx + 1;
    end
    v = -(s(:,rhs_idx)-s(:,lhs_idx))/2/dx;
    
    
    if mod(Ny,2) == 1
        % odd number of grid points
        bot_idx = (Ny-1)/2;
        top_idx = bot_idx + 2;
    else
        % even number of grid points
        bot_idx = Ny/2;
        top_idx = bot_idx + 1;
    end
    u = (s(top_idx,:)-s(bot_idx,:))/2/dy;
    
    xgrid = linspace(0,Lx,Nx);
    ygrid = linspace(0,Ly,Ny);
    
    figure(1)
    hold on
    plot(xgrid,u);
    
    figure(2)
    hold on
    plot(ygrid,v);
    
    if i==1
        [xcoord,ycoord] = meshgrid(xgrid,ygrid);
        figure(3)
        contourf(xcoord,ycoord,w,100);
        
        figure(4)
        contourf(xcoord,ycoord,s);
    end
    
    xmin(i) = xcoord(s==min(s,[],'all'));
    ymin(i) = ycoord(s==min(s,[],'all'));
end

figure(1)
title('Plot of $u$ agaisnt $x$ for $y=0.5,\,L_x=1,\,L_y=1$ on a 161x161 grid');
grid on
grid minor
xlabel("x");
ylabel("u");
legend(filenames(:,2), "Location", "southwest");
saveas(gcf,'Images/u_vs_x.png');

figure(2)
title('Plot of $v$ agaisnt $y$ for $x=0.5,\,L_x=1,\,L_y=1$ on a 161x161 grid');
grid on
grid minor
xlabel("y");
ylabel("v");
legend(filenames(:,2), "Location", "northwest");
saveas(gcf,'Images/v_vs_y.png');

figure(3)
title('Contour plot for vorticity at $Re=100,\,L_x=1,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
saveas(gcf,'Images/w_100.png');

figure(4)
title('Contour plot for streamfunction at $Re=100,\,L_x=1,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
saveas(gcf,'Images/s_100.png');


lx1ly2 = "1_2_161_161_100_data.txt";
lx2ly1 = "2_1_161_161_100_data.txt";

filenames2 = [lx1ly2, lx2ly1];

for i=1:length(filenames2)
    filename = filenames2(i);
    fid = fopen(filename);
    headers = fscanf(fid,"%s \n ",[1,7]);
    headersData = fscanf(fid,"%f, ",[1,7]);
    Lx = headersData(1);
    Ly = headersData(2);
    Nx = headersData(3);
    Ny = headersData(4);
    Re = headersData(7);
    
    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    
    fscanf(fid,"%*s",1);
    w = fscanf(fid,"%f",[Ny,Nx])';
    fscanf(fid,"%*s",1);
    s = fscanf(fid,"%f",[Ny,Nx])';
    
    xgrid = linspace(0,Lx,Nx);
    ygrid = linspace(0,Ly,Ny);

    figure(5+i)
    contourf(xcoord,ycoord,s);
end

figure(6)
title('Contour plot for vorticity at $Re=100,\,L_x=1,\,L_y=2$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
saveas(gcf,'Images/lx1ly2.png');

figure(7)
title('Contour plot for streamfunction at $Re=100,\,L_x=2,\,L_y=1$ on a 161x161 grid');
xlabel("x");
ylabel("y");
colorbar("eastoutside");
saveas(gcf,'Images/lx2ly1.png');
