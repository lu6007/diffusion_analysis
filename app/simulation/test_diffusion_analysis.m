% The unit test function
function test_diffusion_analysis()
cell_name = 'photobleach_cell';
data = sample_diffusion_init_data(cell_name);
data = computer_simulation(data, 'dt', 0.25, 'simulation_time', 40); 
compute_recovery_curve(data); 
pause_str = '\nFunction test_diffusion_analysis(): paused. \nPress any key to close current figures and continue.\n';
fprintf(pause_str); 
pause; 
close all; clear data; 

%
data = sample_diffusion_init_data('layered_diffusion'); 
computer_simulation(data, 'dt', 0.25);
%
data = sample_diffusion_init_data('tensor_diffusion'); 
computer_simulation(data, 'dt', 0.25); 
fprintf(pause_str); 
pause; 
close all; clear data; 

%
cell_name = 'mem17';
data = sample_diffusion_init_data(cell_name);
estimate_frap(data);
statistics_0905_2006;
fprintf(pause_str); 
pause; 
close all; clear data; 

%
cell_name = 'egf_pp1_lyn1'; 
data = sample_diffusion_init_data(cell_name);
correct_fret(data); 

return