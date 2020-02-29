function abercoef = seidelcoef (lens1, lens2, gaussdata)

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
  h = gaussdata.h;
  k = gaussdata.k;
  K = gaussdata.K;
  L = gaussdata.L;
  M = gaussdata.M;
  Md = gaussdata.Md;


  %% ------------- 変数定義 ------------- %%
  N = size(r,2);
  B = zeros(1,N);               % 球面収差係数
  C = zeros(1,N);               % 非点収差係数
  D = zeros(1,N);               % 像面湾曲収差係数
  E = zeros(1,N);               % 歪曲収差係数
  F = zeros(1,N);               % コマ収差係数

  %% ------------- ザイデル収差係数計算 ------------- %%
  for i=1:1:N
    if(i==1)
      N1 = n(1) - n_0;
      N2 = 1/n(1)^2 - 1/n_0^2;
      NS = 1/(n(1)*sd(1)) - 1/(n_0*s(1));
    else
      N1 = n(i) - n(i-1);
      N2 = 1/n(i)^2 - 1/n(i-1)^2;
      NS = 1/(n(i)*sd(i)) - 1/(n(i-1)*s(i));
    end
    B(i) = 1/2*( h(i)^4*(b(i)/r(i)^3)*N1 + h(i)^4*K(i)^2*NS );
    C(i) = 1/2*( h(i)^4*k(i)^2*(b(i)/r(i)^3)*N1 + (1 + h(i)^2*k(i)*K(i))^2*NS );
    D(i) = 1/2*( h(i)^4*k(i)^2*(b(i)/r(i)^3)*N1 + h(i)^2*k(i)*K(i)*(2 + h(i)^2*k(i)*K(i))*NS - K(i)*N2 );
    E(i) = 1/2*( h(i)^4*k(i)^3*(b(i)/r(i)^3)*N1 + k(i)*(1 + h(i)^2*k(i)*K(i))*(2 + h(i)^2*k(i)*K(i))*NS - (1 + h(i)^2*k(i)*K(i))/h(i)^2*N2 );
    F(i) = 1/2*( h(i)^4*k(i)*(b(i)/r(i)^3)*N1 + h(i)^2*K(i)*(1 + h(i)^2*k(i)*K(i))*NS );
  end

  abercoef.B = B;
  abercoef.C = C;
  abercoef.D = D;
  abercoef.E = E;
  abercoef.F = F;

  %% ------------- 色収差計算 ------------- %%
  if(isfield(lens1,'dn'))
    dn = lens1.dn;
    dn_0 = lens2.dn_0;
    ds = zeros(1,N);               % 軸上色収差
    dM = zeros(1,N);               % 色倍率収差
    
    ds = -(1/n(end))*(sd(end)/h(end))^2*h.^2.*K.*(dn./n - [dn_0 dn(1:end-1)]./[n_0 n(1:end-1)]);
    ds_tot = sum(ds);
    dM = (dn_0/n_0 - dn(end)/n(end)) + (ds_tot/(sd(end)-td(end))) - sum(h.^2.*k.*K.*(dn./n - [dn_0 dn(1:end-1)]./[n_0 n(1:end-1)]));

    abercoef.ds = ds;
    abercoef.dM = dM;
  end
  
endfunction
