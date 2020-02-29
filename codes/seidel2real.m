function [DX, DY, Rho, Theta] = seidel2real (lens1, lens2, gaussdata, abercoef, y_0, EPD)

  %% ------------- �����Y�f�[�^���A���o���h�� ------------- %%
  r = lens1.r;                         
  b = lens1.b;           
  d = lens1.d;
  n = lens1.n;
  n_0 = lens2.n_0;
  s_1 = lens2.s_1;
  t_1 = lens2.t_1;

  %% ------------- �K�E�X���w�v�Z�f�[�^���A���o���h�� ------------- %%
  s = gaussdata.s;
  sd = gaussdata.sd;
  t = gaussdata.t;
  td = gaussdata.td;
  M = gaussdata.M;
  Md = gaussdata.Md;

  %% ------------- �����W���f�[�^���A���o���h�� ------------- %%
  B = abercoef.B;
  C = abercoef.C;
  D = abercoef.D;
  E = abercoef.E;
  F = abercoef.F;

  lambda_0 = 1;
  B_tot = sum(B);     % �����W���̊e�ʍ��v
  C_tot = sum(C);
  D_tot = sum(D);
  E_tot = sum(E);
  F_tot = sum(F);

  rho = [0.000001 (0.01:0.01:1)]*(EPD/2);    % �ˏo���ɂ�����a�������� (mm)
  theta = deg2rad([0 90 180 270 (0:15:360)])';
  [Rho, Theta] = meshgrid(rho,theta);
  Dx = B_tot*Rho.^3.*sin(Theta) - 2*F_tot*y_0*Rho.^2.*sin(Theta).*cos(Theta) + D_tot*y_0^2*Rho.*sin(Theta);
  Dy = B_tot*Rho.^3.*cos(Theta) - F_tot*y_0*Rho.^2.*(1+2*cos(Theta).^2) + (2*C_tot+D_tot)*y_0^2*Rho.*cos(Theta) - E_tot*y_0^3;
  DX = (td(end)-sd(end))/(n(end)*Md)*Dx;
  DY = (td(end)-sd(end))/(n(end)*Md)*Dy;

endfunction
