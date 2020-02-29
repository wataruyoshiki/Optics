clear all;
close all;

%% ------------- �p�����[�^��`  ------------- %%
r = [100 500];                % �ȗ����a (mm)
b = [0 0];                    % �񋅖ʌW��
d = [50 100];                 % �ʊԊu (mm)
n = [1.44 1];                 % ���ܗ�

s_1 = -5000;                  % ���ʂ��畨�̖ʂ܂ł̋��� (mm)
t_1 = -0.00001;               % ���ʂ�����˓��܂ł̋��� (mm)
n_0 = 1;                      % ���ʂ��畨�̖ʂɂ�����}���̋��ܗ�
phi = deg2rad((0:2:20));      % ��p(���p) (rad.)
EPD = 100;                    % ���˓����a (mm)

isFigLens = 'y';
isFigAber = 'n';
isRaytrace = 'n';


%% ------------- �ϐ���` ------------- %%
lambda_0 = 1;                       % 1�ȊO�ɐݒ肵�Ȃ����ƁB
N = size(r,2);                      % ���ܖʐ�
N_phi = size(phi,2);                % ��p������
y_0 = n_0*lambda_0*tan(phi);        % ���K�����ꂽ���̖ʂɂ����镨�_����
Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0; % ���̍���


%% ------------- �����Y�f�[�^���o���h�� ------------- %%
lens1.r = r;                       
lens1.b = b;           
lens1.d = d;
lens1.n = n;
lens2.n_0 = n_0;
lens2.s_1 = s_1;
lens2.t_1 = t_1;


%% ------------- �K�E�X���w�v�Z ------------- %%
gaussdata = gauss (lens1,lens2);
s = gaussdata.s;                    % �K�E�X���w�v�Z�f�[�^���A���o���h��
sd = gaussdata.sd;
t = gaussdata.t;
td = gaussdata.td;
h = gaussdata.h;
k = gaussdata.k;
K = gaussdata.K;
L = gaussdata.L;
M = gaussdata.M;
Md = gaussdata.Md;


%% ------------- ���H�}�`�� ------------- %%
if(isFigLens == 'y')
  figure;
  %% ---- ������� ---- %%
  [tmp1, Yobj1, temp2, Yk1, Zk1] = raytrace3d (lens1, lens2, gaussdata, 0, EPD/2, 0);
  lensdraw(lens1, lens2, gaussdata, Yk1, Zk1, EPD, 0, Yobj1,'b-');
  [tmp1, Yobj1, temp2, Yk1, Zk1] = raytrace3d (lens1, lens2, gaussdata, 0, -EPD/2, 0);
  lensdraw(lens1, lens2, gaussdata, Yk1, Zk1, EPD, 0, Yobj1,'b-');
  %% ---- ���O���� ---- %%
  [tmp1, Yobj2, temp2, Yk2, Zk2] = raytrace3d (lens1, lens2, gaussdata, max(y_0), EPD/2, 0);
  lensdraw(lens1, lens2, gaussdata, Yk2, Zk2, EPD, max(y_0), Yobj2,'b--');
  [tmp1, Yobj2, temp2, Yk2, Zk2] = raytrace3d (lens1, lens2, gaussdata, max(y_0), -EPD/2, 0);
  lensdraw(lens1, lens2, gaussdata, Yk2, Zk2, EPD, max(y_0), Yobj2,'b--');
end

%% ------------- �����W���v�Z ------------- %%
[B, C, D, E, F] = seidelcoef(lens1, lens2, gaussdata);


%% ------------- �������ʌv�Z ------------- %%
for j=1:1:N_phi
  [DX(:,:,j), DY(:,:,j), Rho, Theta] = seidel2real (lens1, lens2, gaussdata, B, C, D, E, F, y_0(j), EPD);
end


%% ------------- �������ǐ� ------------- %%
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
end

%% ------------- �c�����v�Z ------------- %%
[DZ_X, DZ_Y, Dist_Y] = lat2lon (DX, DY, gaussdata, Rho, Theta, Y_0);
##[DZ_X_rt, DZ_Y_rt, Dist_Y_rt] = lat2lon (DX_rt, DY_rt, gaussdata, Rho, Theta, Y_0);

%% ------------- �����Ȑ��`�� ------------- %%
  if(isFigAber=='y')
  %% ---- �������Ȑ� ---- %%
  f1 = figure;
  subplot(1,2,1);       % �^���W�F���V����
  plot(Rho(1,:),DY(1,:,1),'b-',Rho(1,:),DY(1,:,fix(N_phi/2+1)),'r-',Rho(1,:),DY(1,:,end),'g-');
  hold on;
  plot(-Rho(3,:),DY(3,:,1),'b-',-Rho(3,:),DY(3,:,fix(N_phi/2+1)),'r-',-Rho(3,:),DY(3,:,end),'g-');
  hold off;
  xlabel('ENP (mm)');
  ylabel('DY (mm)');
  legend([num2str(rad2deg(phi(1))) ' deg.'],[num2str(rad2deg(phi(fix(N_phi/2+1)))) ' deg.'],[num2str(rad2deg(phi(end))) ' deg.']);
  set(gca,'FontSize',16 );

  subplot(1,2,2);       % �T�W�^��
  plot(Rho(2,:),DX(2,:,1),'b-',Rho(2,:),DX(2,:,fix(N_phi/2+1)),'r-',Rho(2,:),DX(2,:,end),'g-');
  xlabel('ENP (mm)');
  ylabel('DX (mm)');
  legend([num2str(rad2deg(phi(1))) ' deg.'],[num2str(rad2deg(phi(fix(N_phi/2+1)))) ' deg.'],[num2str(rad2deg(phi(end))) ' deg.']);
  set(gca,'FontSize',16 );

  %% ---- �c�����Ȑ� ---- %%
  figure;
  subplot(1,3,1);       % ���ʎ���
  plot(DZ_Y(1,:,1),Rho(1,:),'b-');
  ##hold on;
  ##plot(DZ_Y_rt(1,:,1),Rho(1,:),'r-');
  ##hold off;
  xlim([-max(abs(DZ_Y(1,:,1))) max(abs(DZ_Y(1,:,1)))]);
  xlabel('Spherical aber. (mm)');
  ylabel('ENP (mm)');
  set(gca,'FontSize',16 );

  subplot(1,3,2);       % ��_����
  plot(DZ_Y(1,end,:),rad2deg(phi),'b-', DZ_X(2,end,:),rad2deg(phi),'b--');
  ##hold on;
  ##plot(DZ_Y_rt(1,end,:),rad2deg(phi),'r-', DZ_X_rt(2,end,:),rad2deg(phi),'r--');
  ##hold off;
  xlabel('Astigmatism (mm)');
  ylabel('Field angle (deg.)');
  legend('T','S');
  set(gca,'FontSize',16 );

  subplot(1,3,3);       % �c�Ȏ���
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
    subplot(1,2,1);       % �^���W�F���V����
    hold on;
    plot(Rho(1,:),DY_rt(1,:,1),'b--',Rho(1,:),DY_rt(1,:,fix(N_phi/2+1)),'r--',Rho(1,:),DY_rt(1,:,end),'g--');
    plot(-Rho(3,:),DY_rt(3,:,1),'b--',-Rho(3,:),DY_rt(3,:,fix(N_phi/2+1)),'r--',-Rho(3,:),DY_rt(3,:,end),'g--');
    hold off;
    subplot(1,2,2);       % �T�W�^��
    hold on;
    plot(Rho(2,:),DX_rt(2,:,1),'b--',Rho(2,:),DX_rt(2,:,fix(N_phi/2+1)),'r--',Rho(2,:),DX_rt(2,:,end),'g--');
    hold off;
  end

end
