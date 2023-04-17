function x = renormalize(x)
  %{
  Change normalization of energy
  %}
  
  H = hamiltonian(x);
  if ( H > 0 ) 
    %currently a scattering state, don't rescale anything. Hopefully we
    %will find our way back to negative energy.
    return;
  end

  l = sqrt(abs(H));

  x(1:4) = x(1:4)*l^2;
  x(5:8) = x(5:8)/l;
  if numel(x) == 9
    x(9)   = x(9)*l^3;
  end
end