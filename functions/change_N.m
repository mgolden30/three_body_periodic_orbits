function y = change_N(y, N)
  n = (numel(y)-1)/8;
  T = y(8*n+1);

  x0 = reshape( y(1:8*n), [8,n] );
  x  = fft( x0, [], 2 );
  x2 = zeros(8,N);

  if N <= n
    %Only copy positive frequencies
    x2(:,1)     =   x(:,1);
    x2(:,2:N/2) = 2*x(:,2:N/2);
  end

  if N > n
    %Only copy positive frequencies
    x2(:,1)     =   x(:,1);
    x2(:,2:n/2) = 2*x(:,2:n/2);
  end

  xf = real(ifft(x2,[],2))/n*N;
  
  plot( x0(1,1:2:end) );
  hold on
    plot( xf(1,:) );
  hold off
  

  y = [reshape(xf, [8*N,1]); T];
end