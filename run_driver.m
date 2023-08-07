% This file runs a series of simulations with varying parameters and stores
% the data

fd_values    = 1.8 : 0.4 : 3.0;   % Fractal Dimensions
stick_values = 0.1 : 0.2 : 1.0;   % Stickiness

nsec_values = 25;     % Number of size classes

n_runs = length(fd_values) * length(stick_values);

run_data = cell(n_runs, 8);

i_run = 1;

for indx = 1 : length(fd_values)
    
    for jindx = 1 : length(stick_values)
            
        disp([num2str(fd_values(indx)) ' ' num2str(stick_values(jindx))])

        [diam_i, total_flux, total_mass, nspec_i, masspec_i, fluxspec_i, p] = coag_driver(nsec_values, fd_values(indx), stick_values(jindx));
        
        run_data{i_run,1} = fd_values(indx);
        run_data{i_run,2} = stick_values(jindx);
        
        run_data{i_run,3} = diam_i;
        run_data{i_run,4} = total_flux;
        run_data{i_run,5} = total_mass;
        run_data{i_run,6} = nspec_i;
        run_data{i_run,7} = masspec_i;
        run_data{i_run,8} = fluxspec_i;
        

        i_run = i_run + 1;
        

    end


end

save runoutput.mat

