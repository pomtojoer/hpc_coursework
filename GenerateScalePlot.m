clc;
close all;

%% Setting the graph interpreter as latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter','latex');

%% Defining the file names and locations to be plotted
filename = "executionTime.txt";
fid = fopen(filename);
fgetl(fid);

numberofprocs = zeros(1,5);
processtimes = zeros(5,5);
for i=1:5
    numberofprocs(i) = fscanf(fid,"%d process",1);
    processtimes(:,i) = fscanf(fid,"%f seconds",[5,1]);
end

meanprocesstimes = mean(processtimes);
normalisedprocesstimes = processtimes/meanprocesstimes(5);
normalisedmeanprocesstimes = meanprocesstimes/meanprocesstimes(5);
numberofprocs = repmat(numberofprocs,5,1);

hold on
figure(1)
scatter(numberofprocs(:),normalisedprocesstimes(:), 'kx');
plot(numberofprocs(1,:),normalisedmeanprocesstimes, 'r-')
title('Scaling Plot: Normalised process times against number of processes used');
grid on
grid minor
xlabel("Number of processes");
ylabel("Normalised process times (time taken/mean time taken for 1 process)");
legend("Actual Time Taken", "Average Time Taken", "Location", "northeast");
saveas(gcf,'Images/scaleplot.png');
