clear all;
close all;

%% ------------- パラメータ定義  ------------- %%
r = [100 500];                % 曲率半径 (mm)
b = [0 0];                    % 非球面係数
d = [50 100];                 % 面間隔 (mm)
n = [1.44 1];                 % 屈折率

s_1 = -5000;                  % 第一面から物体面までの距離 (mm)
t_1 = -0.00001;               % 第一面から入射瞳までの距離 (mm)
n_0 = 1;                      % 第一面から物体面における媒質の屈折率
phi = deg2rad(0);             % 画角(半角) (rad.)
EPD = 100;                    % 入射瞳直径 (mm)

isLensdraw = 'y';
isRaytrace = 'y';
isSeidel = 'y';
isGraph = 'y';


%% ------------- 変数定義 ------------- %%
lambda_0 = 1;                       % 1以外に設定しないこと。
N = size(r,2);                      % 屈折面数

%% ------------- レンズデータをバンドル ------------- %%
lens1.r = r;                       
lens1.b = b;           
lens1.d = d;
lens1.n = n;
lens2.n_0 = n_0;
lens2.s_1 = s_1;
lens2.t_1 = t_1;


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
y_0 = n_0*lambda_0*tan(phi);        % 正規化された物体面における物点高さ
Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0; % 物体高さ
[tmp, Yk, Zk] = raytrace(lens1, lens2, gaussdata, y_0, EPD/2);
lensdraw(lens1, lens2, gaussdata, Yk, Zk, EPD);


%% ------------- 収差係数計算 ------------- %%
[B, C, D, E, F] = seidelcoef(lens1, lens2, gaussdata);


%% ------------- 実収差量計算 ------------- %%
[DX, DY, Rho, Theta] = seidel2real (lens1, lens2, gaussdata, B, C, D, E, F, y_0, EPD);


%% ------------- 実光線追跡 ------------- %%
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
