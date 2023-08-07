% Script to analyze large output

load run4output.mat

% Want to create plots of the error in the sectional model plotted against
% N0 and v0, with a different plot for different numbers of sections. 

num_d0   = length(d0_values);
num_dz   = length(dz_values);
num_nsec = length(nsec_values);

% Extract the data and plot out fluxes. 

y_error    = zeros(num_d0, num_dz, num_nsec);

for indx = 1 : num_nsec * num_d0 * num_dz
    
    d0   = run_data{indx, 1};
    d0   = find(d0 == d0_values);
    dz   = run_data{indx, 2};
    dz   = find(dz == dz_values);
    
    
        
end

