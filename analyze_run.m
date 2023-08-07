% Script to analyze large output

load run2output.mat

% Want to create plots of the error in the sectional model plotted against
% N0 and v0, with a different plot for different numbers of sections. 

num_N0   = length(N0_values);
num_v0   = length(v0_values);
num_nsec = length(nsec_values);

% Calculate the errors between the sectional and analytical values. 

y_error    = zeros(num_N0, num_v0, num_nsec);

for indx = 1 : num_nsec * num_N0 * num_v0
    
    N0   = run_data{indx, 1};
    N0   = find(N0 == N0_values);
    v0   = run_data{indx, 2};
    v0   = find(v0 == v0_values);
    nsec = run_data{indx, 3};
    nsec = find(nsec == nsec_values);
    
    y      = run_data{indx, 5};
    y_quad = run_data{indx, 6};
     
    y_error(N0, v0, nsec) = max(max(abs((y_quad - y)./y_quad)*100));
    
end

% Now work through each set of (N0, v0) pairs for each nsec value and plot
% the surface of error values

N0_mat = N0_values';
N0_mat = N0_mat(:, ones(1, num_v0));

v0_mat = v0_values(ones(num_N0,1), :); 

for i_sec = 1 : num_nsec
    
    figure(i_sec)
    surf(N0_mat, v0_mat, log10(reshape(y_error(:,:,i_sec), 10, 10)));
    xlabel('N_0')
    ylabel('v_0')
    title(['Number of sections = ' num2str(nsec_values(i_sec))])
    
end
