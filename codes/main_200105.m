clear all;
close all;
pkg load statistics;

 
%% ------------- パラメータ設定  ------------- %%
l_F = 486.1;
l_D = 589.2;
l_C = 656.3;

glass_a = 'NBK-7';
glass_b = 'F2';
n_a = RIcalc(glass_a,l_D);
n_b = RIcalc(glass_b,l_D);

abbe_a = (RIcalc(glass_a,l_D)-1)/(RIcalc(glass_a,l_F)-RIcalc(glass_a,l_C));
abbe_b = (RIcalc(glass_b,l_D)-1)/(RIcalc(glass_b,l_F)-RIcalc(glass_b,l_C));

f = 100;      % 全系焦点距離 (mm)
P = 1/100;    % 全系パワー (1/mm)

%% ------------- パワー配分決定 ------------- %%
P_mat = inv([1/abbe_a 1/abbe_b; 1 1])*[0; P];
P_a = P_mat(1);
P_b = P_mat(2);

lens1.n = [n_a n_b 1];
lens1.dn = [RIcalc(glass_a,l_F)-RIcalc(glass_a,l_C) RIcalc(glass_b,l_F)-RIcalc(glass_b,l_C) 0];
lens1.b = [0 0 0];           
lens1.d = [0 0 0];
lens2.n_0 = 1;
lens2.dn_0 = 0;
lens2.s_1 = -500000;
conf.isFigLens = 'n';
conf.isFigAber = 'n';
conf.isRaytrace = 'n';
theta = deg2rad(0);                   % 画角(半角) (rad.)
EPD = 50;                             % 入射瞳直径 (mm)


r_3 = [(-400:10:-50) (50:10:400)];
r_2 = (P_b/(n_b-1) + 1./r_3).^-1;
r_1 = (P_a/(n_a-1) + 1./r_2).^-1;

t_1 = [-0.00001 -5 -10 -15 -20]';
figure;
plot(r_3,r_1,'bo');
hold on; plot(r_3,r_2,'ro');
ylim([-200 200]);
xlabel('r_3 (mm)'); ylabel('r_1, r_2 (mm)');
legend('r_1','r_2');


%% ------------- 収差係数計算 ------------- %%
N = size(r_3,2);
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
    lens1.r = [r_1(i) r_2(i) r_3(i)];                       
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





