function n = RIcalc (glass, lambda)
%% lambda: îgí∑ (nm)
%% glass: è…çﬁñº (ï∂éöóÒ)

  l = lambda*10^-3;       % îgí∑ (um)
  
  if(strcmp(glass,'NBK-7'))
    n_1 = 1.03961212*l^2/(l^2 - 0.00600069867);
    n_2 = 0.231792344*l^2/(l^2 - 0.0200179144);
    n_3 = 1.01046945*l^2/(l^2 - 103.560653);
    n = sqrt(1 + n_1 + n_2 + n_3);
  elseif(strcmp(glass,'F2'))
    n_1 = 1.34533359*l^2/(l^2 - 0.00997743871);
    n_2 = 0.209073176*l^2/(l^2 - 0.0470450767);
    n_3 = 0.937357162*l^2/(l^2 - 111.886764);
    n = sqrt(1 + n_1 + n_2 + n_3);
  else
    n = 1;
    disp([glass ': No glass data']);
  end

endfunction
