%Search for close passes of random initial guesses

max_iterations = 2^20;
min_R = 0.025;
max_R = 3.5;
close_pass_threshold = 1;
integration_error = 1e-2;
min_h = 1e-9;

search_for_close_passes( max_iterations, min_R, max_R, min_h, close_pass_threshold, integration_error );