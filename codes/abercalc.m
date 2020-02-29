function [aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, phi, EPD, conf)


  %% ------------- レンズデータをアンバンドル ------------- %%
  r = lens1.r;                         
  b = lens1.b;           
  d = lens1.d;
  n = lens1.n;
  n_0 = lens2.n_0;
  s_1 = lens2.s_1;
  t_1 = lens2.t_1;
  isFigLens = conf.isFigLens;
  isFigAber = conf.isFigAber;
  isRaytrace = conf.isRaytrace;
  

  %% ------------- 変数定義 ------------- %%
  lambda_0 = 1;                       % 1以外に設定しないこと。
  N = size(r,2);                      % 屈折面数
  N_phi = size(phi,2);                % 画角条件数
  y_0 = n_0*lambda_0*tan(phi);        % 正規化された物体面における物点高さ
  Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0; % 物体高さ


  %% ------------- ガウス光学計算 ------------- %%
  gaussdata = gauss (lens1,lens2);
  s = gaussdata.s;                    % ガウス光学計算データをアンバンドル
  sd = gaussdata.sd;
  t = gaussdata.t;
  td = gaussdata.td;
  h = gaussdata.h;
  k = gaussdata.k;
  K = gaussdata.K;
  L = gaussdata.L;
  M = gaussdata.M;
  Md = gaussdata.Md;


  %% ------------- 光路図描画 ------------- %%
  if(isFigLens == 'y')
    figure;
    %% ---- 軸上光線 ---- %%
    [tmp1, Yobj1, temp2, Yk1, Zk1] = raytrace3d (lens1, lens2, gaussdata, 0, EPD/2, 0);
    lensdraw(lens1, lens2, gaussdata, Yk1, Zk1, EPD, 0, Yobj1,'b-');
    [tmp1, Yobj1, temp2, Yk1, Zk1] = raytrace3d (lens1, lens2, gaussdata, 0, -EPD/2, 0);
    lensdraw(lens1, lens2, gaussdata, Yk1, Zk1, EPD, 0, Yobj1,'b-');
    %% ---- 軸外光線 ---- %%
    [tmp1, Yobj2, temp2, Yk2, Zk2] = raytrace3d (lens1, lens2, gaussdata, max(y_0), EPD/2, 0);
    lensdraw(lens1, lens2, gaussdata, Yk2, Zk2, EPD, max(y_0), Yobj2,'b--');
    [tmp1, Yobj2, temp2, Yk2, Zk2] = raytrace3d (lens1, lens2, gaussdata, max(y_0), -EPD/2, 0);
    lensdraw(lens1, lens2, gaussdata, Yk2, Zk2, EPD, max(y_0), Yobj2,'b--');
  end

  %% ------------- 収差係数計算 ------------- %%
  abercoef = seidelcoef(lens1, lens2, gaussdata);
##  abercoef.B = B;
##  abercoef.C = C;
##  abercoef.D = D;
##  abercoef.E = E;
##  abercoef.F = F;


  %% ------------- 実収差量計算 ------------- %%
  for j=1:1:N_phi
    [DX(:,:,j), DY(:,:,j), Rho, Theta] = seidel2real (lens1, lens2, gaussdata, abercoef, y_0(j), EPD);
  end
  aberdata.DX = DX;
  aberdata.DY = DY;
  aberdata.Rho = Rho;
  aberdata.Theta = Theta; 


  %% ------------- 実光線追跡 ------------- %%
  if(isRaytrace=='y')
    N_Theta = size(Theta,1);
    N_Rho = size(Rho,2);
    DX_rt = zeros(N_Theta,N_Rho,N_phi);
    DY_rt = zeros(N_Theta,N_Rho,N_phi);
    for j=1:1:N_phi
      for i=1:1:N_Theta
        [X_rt_tmp, Y_rt_tmp] = raytrace3d(lens1, lens2, gaussdata, y_0(j), Rho(1,:), Theta(i,1));
        DX_rt(i,:,j) = X_rt_tmp';
        DY_rt(i,:,j) = Y_rt_tmp' - M*Y_0(j);
      end
    end
    aberdata.DX_rt = DX_rt;
    aberdata.DY_rt = DY_rt;
  end

  %% ------------- 縦収差計算 ------------- %%
  [DZ_X, DZ_Y, Dist_Y] = lat2lon (DX, DY, gaussdata, Rho, Theta, Y_0);
  aberdata.DZ_X = DZ_X;
  aberdata.DZ_Y = DZ_Y;
  aberdata.Dist_Y = Dist_Y;
  ##[DZ_X_rt, DZ_Y_rt, Dist_Y_rt] = lat2lon (DX_rt, DY_rt, gaussdata, Rho, Theta, Y_0);

  
  %% ------------- 収差曲線描画 ------------- %%
    if(isFigAber=='y')
    %% ---- 横収差曲線 ---- %%
    f1 = figure;
    subplot(1,2,1);       % タンジェンシャル
    plot(Rho(1,:),DY(1,:,1),'b-',Rho(1,:),DY(1,:,fix(N_phi/2+1)),'r-',Rho(1,:),DY(1,:,end),'g-');
    hold on;
    plot(-Rho(3,:),DY(3,:,1),'b-',-Rho(3,:),DY(3,:,fix(N_phi/2+1)),'r-',-Rho(3,:),DY(3,:,end),'g-');
    hold off;
    xlabel('ENP (mm)');
    ylabel('DY (mm)');
    legend([num2str(rad2deg(phi(1))) ' deg.'],[num2str(rad2deg(phi(fix(N_phi/2+1)))) ' deg.'],[num2str(rad2deg(phi(end))) ' deg.']);
    set(gca,'FontSize',16 );

    subplot(1,2,2);       % サジタル
    plot(Rho(2,:),DX(2,:,1),'b-',Rho(2,:),DX(2,:,fix(N_phi/2+1)),'r-',Rho(2,:),DX(2,:,end),'g-');
    xlabel('ENP (mm)');
    ylabel('DX (mm)');
    legend([num2str(rad2deg(phi(1))) ' deg.'],[num2str(rad2deg(phi(fix(N_phi/2+1)))) ' deg.'],[num2str(rad2deg(phi(end))) ' deg.']);
    set(gca,'FontSize',16 );

    %% ---- 縦収差曲線 ---- %%
    figure;
    subplot(1,3,1);       % 球面収差
    plot(DZ_Y(1,:,1),Rho(1,:),'b-');
    ##hold on;
    ##plot(DZ_Y_rt(1,:,1),Rho(1,:),'r-');
    ##hold off;
    xlim([-max(abs(DZ_Y(1,:,1))) max(abs(DZ_Y(1,:,1)))]);
    xlabel('Spherical aber. (mm)');
    ylabel('ENP (mm)');
    set(gca,'FontSize',16 );

    subplot(1,3,2);       % 非点収差
    plot(DZ_Y(1,end,:),rad2deg(phi),'b-', DZ_X(2,end,:),rad2deg(phi),'b--');
    ##hold on;
    ##plot(DZ_Y_rt(1,end,:),rad2deg(phi),'r-', DZ_X_rt(2,end,:),rad2deg(phi),'r--');
    ##hold off;
    xlabel('Astigmatism (mm)');
    ylabel('Field angle (deg.)');
    legend('T','S');
    set(gca,'FontSize',16 );

    subplot(1,3,3);       % 歪曲収差
    plot(Dist_Y(1,1,:),rad2deg(phi),'b-');
    ##hold on;
    ##plot(Dist_Y_rt(1,1,:),rad2deg(phi),'r-');
    ##hold off;
    xlim([-max(abs(Dist_Y(1,1,:))) max(abs(Dist_Y(1,1,:)))]);
    xlabel('Distortion (%)');
    ylabel('Field angle (deg.)');
    set(gca,'FontSize',16 );

    if(isRaytrace == 'y')
      figure(f1);
      subplot(1,2,1);       % タンジェンシャル
      hold on;
      plot(Rho(1,:),DY_rt(1,:,1),'b--',Rho(1,:),DY_rt(1,:,fix(N_phi/2+1)),'r--',Rho(1,:),DY_rt(1,:,end),'g--');
      plot(-Rho(3,:),DY_rt(3,:,1),'b--',-Rho(3,:),DY_rt(3,:,fix(N_phi/2+1)),'r--',-Rho(3,:),DY_rt(3,:,end),'g--');
      hold off;
      subplot(1,2,2);       % サジタル
      hold on;
      plot(Rho(2,:),DX_rt(2,:,1),'b--',Rho(2,:),DX_rt(2,:,fix(N_phi/2+1)),'r--',Rho(2,:),DX_rt(2,:,end),'g--');
      hold off;
    end

  end

endfunction
