% This script analyzes the results from run_driver
%
%
% The data is organized into a cell array. The organization of the cell
% array is
%
%  run_data{i,1} = fractal dimension
%  run_data{i_run,2} = stick_values(jindx);
%  run_data{i_run,3} = diam_i;
%  run_data{i_run,4} = total_flux;
%  run_data{i_run,5} = total_mass;
%  run_data{i_run,6} = nspec_i;
%  run_data{i_run,7} = masspec_i;
%  run_data{i_run,8} = fluxspec_i;
%

% Clean the workspace
close all
clear all

% Load the data
load runoutput.mat

n_runs = size(run_data,1);

% For each stickiness value, plot the different fractal dimension cases on
% one plot

n_fractal = 4;
n_stick   = 5;

f1 = figure(1);
for i_frac = 1 : n_fractal
    
    subplot(2,2,i_frac)
    loglog(run_data{(i_frac-1)*n_stick+1,3}, run_data{(i_frac-1)*n_stick+1, 6}(end,:), ...
           run_data{(i_frac-1)*n_stick+2,3}, run_data{(i_frac-1)*n_stick+2, 6}(end,:), ...
           run_data{(i_frac-1)*n_stick+3,3}, run_data{(i_frac-1)*n_stick+3, 6}(end,:), ...
           run_data{(i_frac-1)*n_stick+4,3}, run_data{(i_frac-1)*n_stick+4, 6}(end,:), ...
           run_data{(i_frac-1)*n_stick+5,3}, run_data{(i_frac-1)*n_stick+5, 6}(end,:))
    xlabel('Diameter [cm]')
    ylabel('Number Spectrum [cm^{-4}]')
    legend('\alpha = 0.1', '\alpha = 0.3', '\alpha = 0.5', '\alpha = 0.7', '\alpha = 0.9', 'Location', 'southwest')
    title(['Steady State Spectra: D = ' num2str(run_data{(i_frac-1)*n_stick+1,1})]) 
end
orient(f1,'landscape')
print -dpdf figure1.pdf

f2 = figure(2);
for i_stick = 1 : n_stick
    
    subplot(5,1,i_stick)
    loglog(run_data{i_stick,3}, run_data{i_stick, 6}(end,:), ...
           run_data{(i_stick + n_stick),3}, run_data{(i_stick + n_stick), 6}(end,:), ...
           run_data{(i_stick + 2*n_stick),3}, run_data{(i_stick + 2*n_stick), 6}(end,:), ...
           run_data{(i_stick + 3*n_stick),3}, run_data{(i_stick + 3*n_stick), 6}(end,:))
    xlabel('Diameter [cm]')
    ylabel('Number Spectrum [cm^{-4}]')
    legend('D = 1.8', 'D = 2.2', 'D = 2.6', 'D = 3.0')
    title(['Steady State Spectra: \alpha = ' num2str(run_data{i_stick,2})]) 
end
orient(f2,'tall')
print -dpdf figure2.pdf

f3 = figure(3);
for i_frac = 1 : n_fractal
    time = [0 : 0.01 : 1.0];
    subplot(2,2,i_frac)
    plot(time, run_data{(i_frac-1)*n_stick+1,4}/10, ... 
         time, run_data{(i_frac-1)*n_stick+2,4}/10, ... 
         time, run_data{(i_frac-1)*n_stick+3,4}/10, ... 
         time, run_data{(i_frac-1)*n_stick+4,4}/10, ... 
         time, run_data{(i_frac-1)*n_stick+5,4}/10)
     xlabel('Time [d]')
     ylabel('Total Flux [mg m^{-2} d^{-1}')
     legend('\alpha = 0.1', '\alpha = 0.3', '\alpha = 0.5', '\alpha = 0.7', '\alpha = 0.9', 'Location', 'northwest')
     title(['Total Flux: D = ' num2str(run_data{(i_frac-1)*n_stick+1,1})])
end
orient(f3,'landscape')
print -dpdf figure3.pdf

f4 = figure(4);
for i_stick = 1 : n_stick
    time = [0 : 0.01 : 1.0];
    subplot(n_stick,1,i_stick)
    plot(time, run_data{i_stick,4}/10, ... 
         time, run_data{(i_stick + n_stick),4}/10, ... 
         time, run_data{(i_stick + 2*n_stick),4}/10, ... 
         time, run_data{(i_stick + 3*n_stick),4}/10)
     xlabel('Time [d]')
     ylabel('Total Flux [mg m^{-2} d^{-1}')
     legend('D = 1.8', 'D = 2.2', 'D = 2.6', 'D = 3.0', 'Location', 'northwest')
     title(['Total Flux: \alpha = ' num2str(run_data{i_stick,2})]) 
end
orient(f4,'tall')
print -dpdf figure4.pdf

