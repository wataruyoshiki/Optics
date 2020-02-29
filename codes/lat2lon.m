function [DZ_X, DZ_Y, Dist_Y] = lat2lon (DX, DY, gaussdata, Rho, Theta, Y_0)

  %% ------------- ガウス光学計算データをアンバンドル ------------- %%
  s = gaussdata.s;
  sd = gaussdata.sd;
  t = gaussdata.t;
  td = gaussdata.td;
  M = gaussdata.M;
  Md = gaussdata.Md;
  
  N_Theta = size(Theta,1);
  N_Rho = size(Rho,2);
  N_phi = size(Y_0,2);                

  DZ_X = zeros(2,N_Rho,N_phi);
  DZ_Y = zeros(2,N_Rho,N_phi);
  Dist_Y = zeros(1,1,N_phi);
  
  for j=1:1:N_phi
    X_obj = DX(1:2,:,j) + M*0;
    Y_obj = DY(1:2,:,j) + M*Y_0(j);
    Z_obj = sd(end)-td(end);
    alpha_X = 1/(Z_obj)*(X_obj - Md*Rho(1:2,:).*sin(Theta(1:2,:)));
    alpha_Y = 1/(Z_obj)*(Y_obj - Md*Rho(1:2,:).*cos(Theta(1:2,:)));
    alphad_X = 1/(Z_obj)*(M*0 - Md*0.*sin(Theta(1:2,:)));
    alphad_Y = 1/(Z_obj)*(M*Y_0(j) - Md*0.*cos(Theta(1:2,:)));
##    DZ_X(:,:,j) = -(X_obj)./(alpha_X);
##    DZ_Y(:,:,j) = -(Y_obj)./(alpha_Y);
    DZ_X(:,:,j) = -(X_obj - M*0)./(alpha_X - alphad_X);
    DZ_Y(:,:,j) = -(Y_obj - M*Y_0(j))./(alpha_Y - alphad_Y);
    if(Y_0(j) == 0)
      Dist_Y(1,1,j) = 0; 
    else
      Dist_Y(1,1,j) = (Y_obj(1,1)./(M*Y_0(j)) - 1)*100; 
    end
  end

endfunction
