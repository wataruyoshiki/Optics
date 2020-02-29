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
phi = deg2rad(0);             % ��p(���p) (rad.)
EPD = 100;                    % ���˓����a (mm)

isLensdraw = 'y';
isRaytrace = 'y';
isSeidel = 'y';
isGraph = 'y';


%% ------------- �ϐ���` ------------- %%
lambda_0 = 1;                       % 1�ȊO�ɐݒ肵�Ȃ����ƁB
N = size(r,2);                      % ���ܖʐ�

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
y_0 = n_0*lambda_0*tan(phi);        % ���K�����ꂽ���̖ʂɂ����镨�_����
Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0; % ���̍���
[tmp, Yk, Zk] = raytrace(lens1, lens2, gaussdata, y_0, EPD/2);
lensdraw(lens1, lens2, gaussdata, Yk, Zk, EPD);


%% ------------- �����W���v�Z ------------- %%
[B, C, D, E, F] = seidelcoef(lens1, lens2, gaussdata);


%% ------------- �������ʌv�Z ------------- %%
[DX, DY, Rho, Theta] = seidel2real (lens1, lens2, gaussdata, B, C, D, E, F, y_0, EPD);


%% ------------- �������ǐ� ------------- %%
N_Theta = size(Theta,1);
N_Rho = size(Rho,2);
X_rt = zeros(N_Theta,N_Rho);
Y_rt = zeros(N_Theta,N_Rho);
for i=1:1:N_Theta
  [X_rt_tmp, Y_rt_tmp] = raytrace3d(lens1, lens2, gaussdata, y_0, Rho(1,:), Theta(i,1));
  X_rt(i,:) = X_rt_tmp';
  Y_rt(i,:) = Y_rt_tmp';
end

figure;
plot(Rho(1,:),Y_rt(1,:)-Y_0*M,'r-',Rho(1,:),DY(1,:),'b-');

figure;
plot(X_rt(:,50),Y_rt(:,50)-Y_0*M,'r-',DX(:,50),DY(:,50),'bo');
axis equal;
