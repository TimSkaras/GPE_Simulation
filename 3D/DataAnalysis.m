% close all
clear all

omega = [1.8 1.9 2.0 2.1 2.2];
tmax = [];
kmax = [];
Amin_kmax = [];
hmethod = 1;

% The ordering of the indices corresponds to the omega for each experiment
% e.g., exp 1   omega = 2
%       exp 2   omega = 1.9
%       exp 3   omega = 2.1
%       exp 4   omega = 1.8, T = 400
%       exp 5   omega = 2.2, T = 400
idx = [4 2 1 3 5];

for j = [1]

    filename = strcat('../FaradayExperiment/zsolution_EXP_', string(j), '.txt');
    
if ismember(j, [4 5])
    [k, t, kspace, x, zsolution] = extractKSpace(filename, 9, 400, hmethod);
    
    % Throw out very low modes because they don't correspond to perturbations
    % Also, create column cutoff to ensure we get the first major maximum
    col_cutoff = 1100;
    kspace = kspace(8:end,1:col_cutoff);
    t = t(1:col_cutoff);
    k = k(8:end);

    [row col] = extractMaxIndex(kspace(1:end, 1:end));
    [Amin_row Amin_col] = extractAminMaxIndex(k, t, kspace);
    
elseif ismember(j, [6])
    
    [k, t, kspace, x, zsolution] = extractKSpace(filename, 9, 500, hmethod);
    
    animateZdensity(x, zsolution, 500)
    
    % Throw out very low modes because they don't correspond to perturbations
    % Also, create column cutoff to ensure we get the first major maximum
    col_cutoff = 1100;
    kspace = kspace(8:end,1:end);
    t = t(1:end);
    k = k(8:end);

    [row col] = extractMaxIndex(kspace(1:end, 1:end));
    [Amin_row Amin_col] = extractAminMaxIndex(k, t, kspace);
    
else
    [k, t, kspace, x, zsolution] = extractKSpace(filename, 9, 250, hmethod);
    
%     animateZdensity(x, zsolution, 250)
    
    % Throw out very low modes because they don't correspond to perturbations
    kspace = kspace(8:end,:);
    k = k(8:end);

    [row col] = extractMaxIndex(kspace);
    [Amin_row Amin_col] = extractAminMaxIndex(k, t, kspace);
end


kmax = [kmax k(row)];
tmax = [tmax t(col)];
Amin_kmax = [Amin_kmax k(Amin_row)];

end

% Convert tmax to milliseconds 
tmax = tmax / (2*pi*476) * 1000;
t = t / (2*pi*476) * 1000;

% xdata = x;
% ydata = zsolution(1:end, 800)';
% 
% fitfunc = @(x, xdata) x(1)*real(sqrt(1 - (xdata - 375./2).^2/x(2)^2)).^4;
% 
% x0 = [.01 100];
% % options = optimoptions('lsqcurvefit','FunctionTolerance', 1e-7)
% coeffs = lsqcurvefit(fitfunc, x0, xdata, ydata);



% figure
% hold on
% plot(xdata, ydata)
% plot(xdata, fitfunc(coeffs, xdata))

% figure
% plot(omega, tmax)
% title('Time to Max Amplitude t_{max} vs. Modulation Frequency \omega/\omega_r')
% xlabel('\omega/\omega_r')
% ylabel('t_{max} (ms)')
% axis([1.75 2.25 50 110])

% figure
% hold on
% plot(omega, kmax)
% plot(1.8:.01:2.2, 1.1 * ones(size(1.8:.01:2.2)))
% % plot(omega, Amin_kmax)
% title('Resonant Wave Number k_{max} vs. Modulation Frequency \omega/\omega_r')
% xlabel('\omega/\omega_r')
% ylabel('k_{max}')
% axis( [1.8 2.2 0 2])
% legend('Simulated k_{max}', 'Predicted k_{max}=1.1', 'Amin Criterion')

[K_grid, T_grid] = meshgrid(k, t);
figure
surf(K_grid', T_grid', kspace, 'EdgeColor', 'none'); 
title('Time Evolution of fft')
xlabel('k')
ylabel('time (ms)')

figure
hold on
semilogy(t, kspace(Amin_row + 1,1:end))
semilogy(t, kspace(Amin_row,1:end))
semilogy(t, kspace(Amin_row - 1,1:end))
% semilogy(t, kspace(52,1:end))
xlabel('time')
ylabel('Mode Amplitude')
set(gca, 'YScale', 'log');

% addFourierNoise(x, zsolution)

% 
% figure
% surf(K_grid', T_grid', sqrt(K_grid'.*kspace)); 
% title('Time Evolution of Power Spectrum')
% xlabel('k \Psi_k')
% ylabel('time')


% % Now construct kmax plot using Amin's criterion of domain validity
% power_spec = sqrt(K_grid'.*kspace);
% temp = ones(size(power_spec))./abs(.009 - power_spec);
% 
% figure
% surf(K_grid', T_grid', power_spec); 
% test = abs(power_spec - .009*ones(size(power_spec))) < .001;
% [M I] = max(test);
% [row col] = extractAminMaxIndex(k, t, kspace);


function [k, t, kspace, x, zsolution] = extractKSpace(filename, start_index, time, henrysmethod)
% Takes a filename with zsolution data and outputs fourier transform matrix
% -- a matrix containing the time evolution of each mode
% 
% INPUTS:
%   filename -- path name of file with data
%   start_index -- first column containing data in filename
%   time -- total time solved for in pde solver (for how many units of time
%   was the system evolved
%   henrysmethod -- either 0 or 1 for how to separate perturbations from
%   ground state. 0 if you want the intial condition subtracted from all
%   future time slices or 1 to use henry's paramter fitting scheme
% 
% The kspace matrix has a mode at each index and the k vector tells you the
% wave number associated with the index. t is vector giving the time at
% each row

fileID = fopen(filename,'r');
% fileID = fopen('../Results/zsolution1.txt','r');
formatSpec = '%f';

A = fscanf(fileID,formatSpec);
zsolution = listToMatrix(A);

fileID = fopen('../FaradayExperiment/gs-zdensity.txt','r');
% fileID = fopen('zdensity.temp','r');
formatSpec = '%f %f\n';
sizeGS = [2 Inf];
groundstate = fscanf(fileID, formatSpec, sizeGS);
groundstate = groundstate(2,:)';

% Extract constants from file
space_points = A(3);
time_points = size(zsolution, 2);
z_length = A(start_index - 1);
HZ = z_length/space_points;

Fs = space_points/z_length;            % Sampling frequency                    
T = HZ;                 % Sampling period       
L = space_points;       % Length of signal
x = (0:L-1)*T;          % Time vector
X = zeros(size(zsolution));

% Isolate perturbations from the groundstate
xdata = x;
fitfunc = @(x, xdata) x(1)*real(sqrt(1 - (xdata - 375./2).^2/x(2)^2)).^4;
x0 = [.01 100];
opts = optimset('Display','off');
lb = -ones(size(x0)) * inf ;
ub = ones(size(x0)) * inf;

for j=1:size(zsolution, 2)
    if henrysmethod
        ydata = zsolution(1:end, j)';
        coeffs = lsqcurvefit(fitfunc, x0, xdata, ydata,lb, ub, opts);
        
        X(1:end,j) = zsolution(1:end,j) - fitfunc(coeffs, xdata)';
    else
        X(1:end,j) = zsolution(1:end,j) - groundstate;
    end
end

% figure
% plot(x(1:end),X)
% title('Initial Ground State')
% xlabel('x (length units)')
% ylabel('X(x)')

% Perform fft on data
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
% P1(2:end-1, :) = 2*P1(2:end-1, :);

kspace = P1;
k = 2*pi*Fs*(0:(L/2))/L;
t = time * (0 : time_points - 1)/time_points;


end

function mat = listToMatrix(A)

% We use the third value in the list as the indexing variable
index = A(3);
starting_index = 9; % This is the index of the first actual data entry (not just meta-data)
mat = [];
% Use the index to create a matrix with index number of rows in each column
% Discard the first six values becuase they are not solution values
% The first three values are Nx, Ny, Nz resp. and the latter three are
% Lx, Ly, and Lz respectively
for j = starting_index:size(A, 1)
    mat(mod(j - starting_index, index) + 1 , floor((j - starting_index) / index) + 1 ) = A(j);
end


end

function [row col] = extractMaxIndex(kspace)
% This function takes the kspace matrix and outputs the index of the entry
% with the maximum value

% Remember that each row corresponds to a different wave number
[M I] = max(kspace');
[M2 I2] = max(M);
row = I2;
col = I(I2);


end

function [row col] = extractAminMaxIndex(k, t, kspace)

[K_grid, T_grid] = meshgrid(k, t);

power_spec = sqrt(K_grid'.*kspace);
bool_mat = abs(power_spec - .009*ones(size(power_spec))) < .001;
[M1 I1] = max(bool_mat);
col = find(M1);
col = col(1);
row = I1(col);

end

function animateZdensity(x, zsolution, tmax)
% Takes position vector x, matrix zsolution which is the time evolution of
% the linear z-density, and a number tmax which is the length 
% of simulated time and outputs an animation of it

time_points = size(zsolution, 2);
times = 0:tmax/(time_points - 1):tmax;
zmax = max(max(zsolution));

opts = optimset('Display','off');
xdata = x;    
fitfunc = @(c, xdata) c(1)*real(sqrt(1 - (xdata - 375./2).^2/c(2)^2)).^4;
x0 = [.01 100];
lb = -ones(size(x0)) * inf ;
ub = ones(size(x0)) * inf;

for k = 1 : length(times)
    
    ydata = zsolution(1:end, k)';

    coeffs = lsqcurvefit(fitfunc, x0, xdata, ydata, lb, ub, opts);
    
    y = zsolution(1:end, k);
    plot(x, y, x, fitfunc(coeffs, xdata));
    plot(x, y)
    axis([x(1) x(end) 0 zmax])
    title( sprintf('Time Evolution of Linear Z-Density: t = %.1f', times(k)) );
    xlabel('Position')
    ylabel('\psi^2')
    pause( 1/2^6 );
end

end

function addFourierNoise(x, zsolution)
% This function takes the zsolution, adds noise to the initial condition,
% and then returns that initial condition

z_initial = zsolution(1:end,1);
L = 375;
N = size(x,2);         % Number of sample points
Fs = N/L;       % Sampling freq. = samples/(spatial length or period of measurement)   

close all
% plot(x, z_initial) % Before adding noise

% Add noise to wave numbers with whole number of wavelengths
% Thus, k = n *pi / L, exclude n <= 4 because those k are too small
% We can let n go up to the number of sampling points n_max = 560
% z_initial = zeros(size(z_initial));
for n = 15:(N-20)
    amplitude = normrnd(0, 5*10^-4);
    z_initial = z_initial .*(1 + amplitude*exp(i*2*pi*n/N * (0:N-1)' ));
end

figure
plot(x, abs(z_initial'))

z_initial = z_initial - zsolution(1:end,1);

% Perform fft on data
Y = fft(z_initial);
P2 = abs(Y/N);
P1 = P2(1:N/2+1,:);
% P1(2:end-1, :) = 2*P1(2:end-1, :);

kspace = P1;
k = 2*pi*Fs*(0:(N/2))/N;

figure
plot(k, kspace')

end