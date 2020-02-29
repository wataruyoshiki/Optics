clear all;
close all;

%% ------------- 光学系仕様  ------------- %%
theta = deg2rad((0:10:10));     % 画角(半角) (rad.)
EPD = 50;                     % 入射瞳直径 (mm)
f = EPD*3.8;                  % 全系焦点距離 (mm)


%% ------------- パワー配置検討  ------------- %%
phi = 1/f;
phi_a = (1.05:0.05:2)*phi;
phi_b = - phi_a;
d_a = phi./phi_a.^2;
L = d_a.*(1 - phi_a./phi) + 1/phi;

figure; plot(phi_a/phi,L,'b-', phi_a/phi,d_a,'r-');
xlabel('\phi_a / \phi (-)');
ylabel('L, d_a (mm)');
legend('L','d_a');
set(gca,'FontSize',16 );

phi_a = 1.4*phi;
phi_b = - phi_a;
d_a = phi/(phi_a^2);
d_b = 1/phi*(1-phi_a*d_a);

 
%% ------------- 面形状検討  ------------- %%
n_NBK7 = 1.5007;
lens1.n = [n_NBK7 1 n_NBK7 1];
lens1.b = [0 0 0 0];           
lens1.d = [0 d_a 0 d_b];
lens2.n_0 = 1;
lens2.s_1 = -50000;
lens2.t_1 = -0.00001;
conf.isFigLens = 'n';
conf.isFigAber = 'n';
conf.isRaytrace = 'n';

r_a2 = [(-1000:50:-50) (10:50:1000)];
r_b2 = [(-1000:50:-50) (2:2:30)];
r_a1 = (phi_a/(lens1.n(1)-1) + 1./r_a2).^-1;
r_b1 = (phi_b/(lens1.n(3)-1) + 1./r_b2).^-1;

figure; plot(r_a2, r_a1, 'bo');
xlabel('r_{a2} (mm)');
ylabel('r_{a1} (mm)');
set(gca,'FontSize',16 );

figure; plot(r_b2, r_b1, 'bo');
xlabel('r_{b2} (mm)');
ylabel('r_{b1} (mm)');
set(gca,'FontSize',16 );

N_a = size(r_a2,2);
N_b = size(r_b2,2);
DR1 = zeros(N_a,N_b);
DR2 = zeros(N_a,N_b);
for i_a=1:1:N_a
  disp([num2str(i_a) ' / ' num2str(N_a)]);
  for i_b=1:1:N_b
    lens1.r = [r_a1(i_a) r_a2(i_a) r_b1(i_b) r_b2(i_b)];                       
    [aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, theta, EPD, conf);
    DX_1 = aberdata.DX(5:end,:,1);
    DY_1 = aberdata.DY(5:end,:,1);
    DX_2 = aberdata.DX(5:end,:,end);
    DY_2 = aberdata.DY(5:end,:,end);
    
    [N_X, N_Y] = size(DX_1(:,:,1));
    DX_1_mu = mean(mean(DX_1));
    DX_2_mu = mean(mean(DX_2));
    DY_1_mu = mean(mean(DY_1));
    DY_2_mu = mean(mean(DY_2));
    
    DX_1_sigma = sqrt(mean(mean((DX_1 - DX_1_mu).^2)));
    DX_2_sigma = sqrt(mean(mean((DX_2 - DX_2_mu).^2)));
    DY_1_sigma = sqrt(mean(mean((DY_1 - DY_1_mu).^2)));
    DY_2_sigma = sqrt(mean(mean((DY_2 - DY_2_mu).^2)));
    
    DR1(i_a,i_b) = sqrt(DX_1_sigma^2 + DY_1_sigma^2);
    DR2(i_a,i_b) = sqrt(DX_2_sigma^2 + DY_2_sigma^2);
    
  end
end


[val1, i_a_min1] = min(min(DR1,[],2),[],1);
[val2, i_b_min1] = min(DR1(i_a_min1,:),[],2);
[val3, i_a_min2] = min(min(DR2,[],2),[],1);
[val4, i_b_min2] = min(DR2(i_a_min2,:),[],2);
disp(['r_a1 ' 'r_a2 ' 'r_b1 ' 'r_b2']);
disp([r_a1(i_a_min1) r_a2(i_a_min1) r_b1(i_b_min1) r_b2(i_b_min1)]);
disp([r_a1(i_a_min2) r_a2(i_a_min2) r_b1(i_b_min2) r_b2(i_b_min2)]);
val1
val3

##r = [100 300];                % 曲率半径 (mm)
##b = [0 0];                    % 非球面係数
##d = [80 100];                 % 面間隔 (mm)
##n = [1.44 1];                 % 屈折率
##
##s_1 = -5000;                  % 第一面から物体面までの距離 (mm)
##t_1 = 50;               % 第一面から入射瞳までの距離 (mm)
##n_0 = 1;                      % 第一面から物体面における媒質の屈折率
##
##isFigLens = 'y';
##isFigAber = 'y';
##isRaytrace = 'y';
##
##
##%% ------------- データをバンドル ------------- %%
##lens1.r = r;                       
##lens1.b = b;           
##lens1.d = d;
##lens1.n = n;
##lens2.n_0 = n_0;
##lens2.s_1 = s_1;
##lens2.t_1 = t_1;
##conf.isFigLens = isFigLens;
##conf.isFigAber = isFigAber;
##conf.isRaytrace = isRaytrace;
##
##[aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, phi, EPD, conf);