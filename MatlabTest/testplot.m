filename = 'data.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

subplot(221), contour(((w))), axis('square');            % Vorticity
subplot(222), contour(((s))), axis('square');