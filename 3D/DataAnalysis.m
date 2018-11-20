close all
clear all

omega = [1.8 1.9 2.0 2.1 2.2];
tmax = [];
kmax = [];
Amin_kmax = [];

% The ordering of the indices corresponds to the omega for each experiment
% e.g., exp 1   omega = 2
%       exp 2   omega = 1.9
%       exp 3   omega = 2.1
%       exp 4   omega = 1.8, T = 400
%       exp 5   omega = 2.2, T = 400
idx = [4 2 1 3 5];

for i = [6]

    filename = strcat('../FaradayExperiment/zsolution_EXP_', string(i), '.txt');
    
if ismember(i, [4 5])
    [k, t, kspace] = extractKSpace(filename, 9, 400);
    
    % Throw out very low modes because they don't correspond to perturbations
    % Also, create column cutoff to ensure we get the first major maximum
    col_cutoff = 1100;
    kspace = kspace(8:end,1:col_cutoff);
    t = t(1:col_cutoff);
    k = k(8:end);

    [row col] = extractMaxIndex(kspace(1:end, 1:end));
    [Amin_row Amin_col] = extractAminMaxIndex(k, t, kspace);
    
elseif ismember(i, [6])
    
    [k, t, kspace] = extractKSpace(filename, 9, 500);
    
    % Throw out very low modes because they don't correspond to perturbations
    % Also, create column cutoff to ensure we get the first major maximum
    col_cutoff = 1100;
    kspace = kspace(8:end,1:end);
    t = t(1:end);
    k = k(8:end);

    [row col] = extractMaxIndex(kspace(1:end, 1:end));
    [Amin_row Amin_col] = extractAminMaxIndex(k, t, kspace);
    
else
    [k, t, kspace] = extractKSpace(filename, 9, 250);
    
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
% 
% figure
% semilogy(t, kspace(54,1:end))
% xlabel('time')
% ylabel('Mode Amplitude')
% set(gca, 'YScale', 'log');

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


function [k, t, kspace] = extractKSpace(filename, start_index, time)
% Takes a filename with zsolution data and outputs fourier transform matrix
% -- a matrix containing the time evolution of each mode
% 
% INPUTS:
%   filename -- path name of file with data
%   start_index -- first column containing data in filename
%   time -- total time solved for in pde solver (for how many units of time
%   was the system evolved
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

% Isolate perturbations from the groundstate
for i=1:size(zsolution, 2)
    zsolution(1:end,i) = zsolution(1:end,i) - groundstate;
end

X = zsolution;

% figure
% plot(x(1:end),X)
% title('Initial Ground State')
% xlabel('x (length units)')
% ylabel('X(x)')

% Perform fft on data
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1, :) = 2*P1(2:end-1, :);

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
for i = starting_index:size(A, 1)
    mat(mod(i - starting_index, index) + 1 , floor((i - starting_index) / index) + 1 ) = A(i);
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