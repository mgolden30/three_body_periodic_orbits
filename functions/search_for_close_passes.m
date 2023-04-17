function search_for_close_passes( max_iterations, min_R, max_R, min_h, close_pass_threshold, integration_error )
rng(1); %reproducibility
infs = 0; %number of times integration fails
a    = 1; %candidate label
for i = 1:max_iterations
  i

  x = 2*(2*rand(8,1) - 1);
  while hamiltonian(x) > 0
    x = 2*(2*rand(8,1) - 1);
  end
  x = renormalize(x); %fix H=-1
  T = 5 + 10*rand(1); %pick a random period in [5,15]

  %carefully evolve forward in time with RK23
  xf = rk23(x, T, integration_error, min_R, max_R, min_h);
  
  if(xf == inf)
    infs = infs + 1;
  end

  if( norm(x - xf) < close_pass_threshold )
    x = [x;T];
    save("candidates/" + a, "x");
    a = a+1;
  end
end
fprintf("%d infs found from integration\n", infs);
end


function x = rk23(x, T, error, min_R, max_R, min_h)
  t = 0; %time integrated
  h = 0.1; %guess of timestep

  k1 = h*f(x);
  while( t < T )
    if( t+h > T )
      h = T-t;
    end

    r1 = x(1:2);
    r2 = x(3:4);
    r12 = r1-r2;
    r13 = 2*r1+r2;
    r23 = 2*r2+r1;

    rs = vecnorm([r12,r23,r13]);
    if (min(rs) < min_R) || (max(rs) > max_R)
      x = inf;
      return;
    end

    k2 = h*f(x+k1/2);
    k3 = h*f(x+k2*3/4);

    x_3rd = x + (2*k1 + 3*k2 + 4*k3)/9; %Ralston third order guess

    k4 = h*f(x_3rd);

    x_2nd = x + (1/4*k1 + 3/10*k2 + 2/5*k3 + 1/20*k4);

    if( norm(x_2nd - x_3rd) > error )
      h = h/2;
      if h < min_h
        x = inf;
        return;
      end
      continue;
    end

    t = t + h;
    h = h*1.01; %increase step size
    x = x_3rd;
    k1= k4;
  end
end