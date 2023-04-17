addpath("functions")

N = 1024;

y = zeros(8*N+1,1);
for i = 1:194
  load( "candidates/" + i + ".mat" );

  y(1:8) = x(1:8);
  y(end) = x(9);

  y = generate_timeseries(y);

  plot_state(y,0);
  drawnow

  saveas(gcf, "candidates/figures/" + i + ".png");              
end