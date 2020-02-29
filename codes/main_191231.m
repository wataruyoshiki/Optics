clear all;
close all;
pkg load statistics;

 
%% ------------- パラメータ設定  ------------- %%
n_NBK7 = 1.5007;
lens1.n = [n_NBK7 1];
lens1.b = [0 0];           
lens1.d = [0 0];
lens2.n_0 = 1;
lens2.s_1 = -50000;
conf.isFigLens = 'n';
conf.isFigAber = 'n';
conf.isRaytrace = 'n';
theta = deg2rad(0);           % 画角(半角) (rad.)
EPD = 50;                     % 入射瞳直径 (mm)


f_1 = 100;
n_1 = n_NBK7;
r_2 = [(-400:10:-50) (50:10:400)];
r_1 = (1/((n_1-1)*f_1) + 1./r_2).^-1;
t_1 = [-0.00001 -5 -10 -15 -20]';
figure;
plot(r_2,r_1,'bo'); ylim([-200 200]);
xlabel('r_2 (mm)'); ylabel('r_1 (mm)');


%% ------------- 収差係数計算 ------------- %%
N = size(r_2,2);
M = size(t_1,1);
B = zeros(M,N);
C = zeros(M,N);
D = zeros(M,N);
E = zeros(M,N);
F = zeros(M,N); 

for j=1:1:M
  disp([num2str(j) ' / ' num2str(M)]);
  lens2.t_1 = t_1(j);
  for i=1:1:N
    lens1.r = [r_1(i) r_2(i)];                       
    [aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, theta, EPD, conf);
    B(j,i) = sum(abercoef.B);
    C(j,i) = sum(abercoef.C);
    D(j,i) = sum(abercoef.D);
    E(j,i) = sum(abercoef.E);
    F(j,i) = sum(abercoef.F);
  end
end



##%% ------------- グラフ描画 ------------- %%
##figure;
##subplot(5,1,1); plot(r_2,B,'bo');
##xlabel('r_2 (mm)'); title('Spherical aber.');
##subplot(5,1,2); plot(r_2,F,'bo');
##xlabel('r_2 (mm)'); title('Coma aber.');
##subplot(5,1,3); plot(r_2,C,'bo');
##xlabel('r_2 (mm)'); title('Astigmatism');
##subplot(5,1,4); plot(r_2,D,'bo');
##xlabel('r_2 (mm)'); title('Field curvature');
##subplot(5,1,5); plot(r_2,E,'bo');
##xlabel('r_2 (mm)'); title('Distortion');





