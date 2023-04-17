function x = generate_timeseries( x )
  N = (numel(x)-1)/8;
  dt= x(end)/N;
  for i = 2:N
    x( 8*(i-1) + (1:8) ) = verlet( x(8*(i-2) + (1:8) ), dt );
  end
end


function x = verlet(x, h)
  x(1:4) = x(1:4) + h/2 * x(5:8);
  v = f(x);
  x(5:8) = x(5:8) + h*v(5:8);
  x(1:4) = x(1:4) + h/2 * x(5:8);
end