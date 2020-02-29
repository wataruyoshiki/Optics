function [X_rt, Y_rt, X, Y, Z] = raytrace3d (lens1, lens2, gaussdata, y_0, rho, theta)

  %% ------------- レンズデータをアンバンドル ------------- %%
  r = lens1.r;                         
  b = lens1.b;           
  d = lens1.d;
  n = lens1.n;
  n_0 = lens2.n_0;
  s_1 = lens2.s_1;
  t_1 = lens2.t_1;
  
  %% ------------- ガウス光学計算データをアンバンドル ------------- %%
  s = gaussdata.s;
  sd = gaussdata.sd;
  t = gaussdata.t;
  td = gaussdata.td;

  lambda_0 = 1;
  eta_0 = rho';
  Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0;                 % 物体高さ
  Xd_0 = eta_0*sin(theta);
  Yd_0 = eta_0*cos(theta);
  Zd_0 = t_1;
  Xc_0 = 0;
  Yc_0 = Y_0;
  Zc_0 = s_1-t_1;
  norm = sqrt((Xd_0-Xc_0).^2 + (Yd_0-Yc_0).^2 + (Zd_0-Zc_0).^2);
  Ld_0 = (Xd_0-Xc_0)./norm;
  Md_0 = (Yd_0-Yc_0)./norm;
  Nd_0 = (Zd_0-Zc_0)./norm;
  
  N2 = size(r,2);           % レンズ面数
  N1 = size(rho',1);        % 瞳径の条件数
  X = zeros(N1,N2);         % 面iへの入射光線が面iと交わる点のX座標(面i座標系)
  Y = zeros(N1,N2);         % 面iへの入射光線が面iと交わる点のY座標(面i座標系)
  Z = zeros(N1,N2);         % 面iへの入射光線が面iと交わる点のZ座標(面i座標系)
  Xc = zeros(N1,N2);        % 面iへの入射光線が面iと交わる点のX座標(面i+1座標系)
  Yc = zeros(N1,N2);        % 面iへの入射光線が面iと交わる点のY座標(面i+1座標系)
  Zc = zeros(N1,N2);        % 面iへの入射光線が面iと交わる点のZ座標(面i+1座標系)
  L = zeros(N1,N2);         % 面iへの入射光線のX方向余弦
  M = zeros(N1,N2);         % 面iへの入射光線のY方向余弦
  N = zeros(N1,N2);         % 面iへの入射光線のZ方向余弦
  Ln = zeros(N1,N2);        % 面iへの入射点における面法線ベクトルのX方向余弦
  Mn = zeros(N1,N2);        % 面iへの入射点における面法線ベクトルのY方向余弦
  Nn = zeros(N1,N2);        % 面iへの入射点における面法線ベクトルのZ方向余弦
  Ld = zeros(N1,N2);        % 面iからの屈折光線のX方向余弦
  Md = zeros(N1,N2);        % 面iからの屈折光線のY方向余弦
  Nd = zeros(N1,N2);        % 面iからの屈折光線のZ方向余弦
  cosI = zeros(N1,N2);      % 面iへの入射角
  cosId = zeros(N1,N2);     % 面iからの屈折角

  
  %% ------------- 3D実光線追跡 ------------- %%
  for i=1:1:N2
    if(i==1)

      L(:,i) = Ld_0;
      M(:,i) = Md_0;
      N(:,i) = Nd_0;
      F = N(:,i) - 1/r(i).*(L(:,i).*Xc_0 + M(:,i).*Yc_0 + N(:,i).*Zc_0);
      G = 1/r(i)*(Xc_0.^2 + Yc_0.^2 + Zc_0.^2) - 2*Zc_0;
      D = r(i)*(F - sqrt(F.^2 - 1/r(i).*G));
      cosI(:,i) = sqrt(F.^2 - 1/r(i).*G);
      X(:,i) = Xc_0 + L(:,i).*D;
      Y(:,i) = Yc_0 + M(:,i).*D;
      Z(:,i) = Zc_0 + N(:,i).*D;

      cosId(:,i) = 1/n(i)*sqrt(n(i)^2 - n_0^2 + n_0^2*cosI(:,i).^2);
      sigma = 1/r(i)*(n(i)*cosId(:,i) - n_0*cosI(:,i));
      Ld(:,i) = 1/n(i)*(n_0*L(:,i) - sigma.*X(:,i));
      Md(:,i) = 1/n(i)*(n_0*M(:,i) - sigma.*Y(:,i));
      Nd(:,i) = 1/n(i)*(n_0*N(:,i) - sigma.*(Z(:,i)-r(i)));

    else

      L(:,i) = Ld(:,i-1);
      M(:,i) = Md(:,i-1);
      N(:,i) = Nd(:,i-1);
      Xc(:,i-1) = X(:,i-1);
      Yc(:,i-1) = Y(:,i-1);
      Zc(:,i-1) = Z(:,i-1) - d(i-1);
      F = N(:,i) - 1/r(i).*(L(:,i).*Xc(:,i-1) + M(:,i).*Yc(:,i-1) + N(:,i).*Zc(:,i-1));
      G = 1/r(i)*(Xc(:,i-1).^2 + Yc(:,i-1).^2 + Zc(:,i-1).^2) - 2*Zc(:,i-1);
      D = r(i)*(F - sqrt(F.^2 - 1/r(i).*G));
      cosI(:,i) = sqrt(F.^2 - 1/r(i).*G);
      X(:,i) = Xc(:,i-1) + L(:,i).*D;
      Y(:,i) = Yc(:,i-1) + M(:,i).*D;
      Z(:,i) = Zc(:,i-1) + N(:,i).*D;

      cosId(:,i) = 1/n(i)*sqrt(n(i)^2 - n(i-1)^2 + n(i-1)^2*cosI(:,i).^2);
      sigma = 1/r(i)*(n(i)*cosId(:,i) - n(i-1)*cosI(:,i));
      Ld(:,i) = 1/n(i)*(n(i-1)*L(:,i) - sigma.*X(:,i));
      Md(:,i) = 1/n(i)*(n(i-1)*M(:,i) - sigma.*Y(:,i));
      Nd(:,i) = 1/n(i)*(n(i-1)*N(:,i) - sigma.*(Z(:,i)-r(i)));

    end
  
    
  end
  
  %% ------------- 像面における光線到達位置 ------------- %%
  Xc_end = X(:,end);
  Yc_end = Y(:,end);
  Zc_end = Z(:,end)-sd(end);
  Z_obj = 0;
  D = (Z_obj - Zc_end)./Nd(:,end);
  X_obj = Xc_end + D.*Ld(:,end);
  Y_obj = Yc_end + D.*Md(:,end);
  
  X_rt = X_obj;
  Y_rt = Y_obj;

endfunction
