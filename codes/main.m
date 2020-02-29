clear all;
close all;

%% ------------- パラメータ定義  ------------- %%
r = [73.5017 -900.0000 9.0675 8.0000];                % 曲率半径 (mm)
b = [0 0 0 0];                    % 非球面係数
d = [0 96.939 0 54.286];                 % 面間隔 (mm)
n = [1.5007 1 1.5007 1];                 % 屈折率

s_1 = -50000;                  % 第一面から物体面までの距離 (mm)
t_1 = -0.0001;               % 第一面から入射瞳までの距離 (mm)
n_0 = 1;                      % 第一面から物体面における媒質の屈折率
phi = deg2rad((0:1:10));      % 画角(半角) (rad.)
EPD = 50;                    % 入射瞳直径 (mm)

isFigLens = 'y';
isFigAber = 'y';
isRaytrace = 'y';


%% ------------- データをバンドル ------------- %%
lens1.r = r;                       
lens1.b = b;           
lens1.d = d;
lens1.n = n;
lens2.n_0 = n_0;
lens2.s_1 = s_1;
lens2.t_1 = t_1;
conf.isFigLens = isFigLens;
conf.isFigAber = isFigAber;
conf.isRaytrace = isRaytrace;

[aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, phi, EPD, conf);