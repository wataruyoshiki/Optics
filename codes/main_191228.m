clear all;
close all;
pkg load statistics;

%% ------------- 光学系仕様  ------------- %%
theta = deg2rad((0:10:20));     % 画角(半角) (rad.)
EPD = 50;                     % 入射瞳直径 (mm)


##%% ------------- パワー配置検討  ------------- %%
##phi = 1/f;
##phi_a = (1.05:0.05:2)*phi;
##phi_b = - phi_a;
##d_a = phi./phi_a.^2;
##L = d_a.*(1 - phi_a./phi) + 1/phi;
##
##figure; plot(phi_a/phi,L,'b-', phi_a/phi,d_a,'r-');
##xlabel('\phi_a / \phi (-)');
##ylabel('L, d_a (mm)');
##legend('L','d_a');
##set(gca,'FontSize',16 );
##
##phi_a = 1.4*phi;
##phi_b = - phi_a;
##d_a = phi/(phi_a^2);
##d_b = 1/phi*(1-phi_a*d_a);

 
%% ------------- パラメータ設定  ------------- %%
n_NBK7 = 1.5007;
lens1.n = [n_NBK7 1];
lens1.b = [0 0];           
lens1.d = [0 0];
lens2.n_0 = 1;
lens1.r = [100.28 -100];                       

lens2.s_1 = -50000;
lens2.t_1 = -0.00001;
conf.isFigLens = 'y';
conf.isFigAber = 'n';
conf.isRaytrace = 'y';

##r_a2 = [(-1000:50:-50) (10:50:1000)];
##r_b2 = [(-1000:50:-50) (2:2:30)];
##r_a1 = (phi_a/(lens1.n(1)-1) + 1./r_a2).^-1;
##r_b1 = (phi_b/(lens1.n(3)-1) + 1./r_b2).^-1;
##
##figure; plot(r_a2, r_a1, 'bo');
##xlabel('r_{a2} (mm)');
##ylabel('r_{a1} (mm)');
##set(gca,'FontSize',16 );
##
##figure; plot(r_b2, r_b1, 'bo');
##xlabel('r_{b2} (mm)');
##ylabel('r_{b1} (mm)');
##set(gca,'FontSize',16 );

[aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, theta, EPD, conf);
DX_1 = aberdata.DX(5:end,:,1);
DY_1 = aberdata.DY(5:end,:,1);
DX_2 = aberdata.DX(5:end,:,2);
DY_2 = aberdata.DY(5:end,:,2);
DX_3 = aberdata.DX(5:end,:,3);
DY_3 = aberdata.DY(5:end,:,3);
DX_rt_1 = aberdata.DX_rt(5:end,:,1);
DY_rt_1 = aberdata.DY_rt(5:end,:,1);
DX_rt_2 = aberdata.DX_rt(5:end,:,2);
DY_rt_2 = aberdata.DY_rt(5:end,:,2);
DX_rt_3 = aberdata.DX_rt(5:end,:,3);
DY_rt_3 = aberdata.DY_rt(5:end,:,3);
Rho = aberdata.Rho;
Theta = aberdata.Theta;
set(gca,'FontSize',16 );


f1=figure; % DX vs DY
for i=1:1:10
  plot(DX_rt_1(:,10*i+1),DY_rt_1(:,10*i+1),'bo');
  hold on;
  plot(DX_1(:,10*i+1),DY_1(:,10*i+1),'r-');
end
axis equal;
xlabel('\Delta X (mm)');
ylabel('\Delta Y (mm)');
set(gca,'FontSize',16 );
axis equal;

f2=figure; % DX vs DY
for i=1:1:10
  plot(DX_rt_2(:,10*i+1),DY_rt_2(:,10*i+1),'bo');
  hold on;
  plot(DX_2(:,10*i+1),DY_2(:,10*i+1),'r-');
end
axis equal;
xlabel('\Delta X (mm)');
ylabel('\Delta Y (mm)');
set(gca,'FontSize',16 );
hold off;

f3=figure; % DX vs DY
for i=1:1:10
  plot(DX_rt_3(:,10*i+1),DY_rt_3(:,10*i+1),'bo');
  hold on;
  plot(DX_3(:,10*i+1),DY_3(:,10*i+1),'r-');
end
axis equal;
xlabel('\Delta X (mm)');
ylabel('\Delta Y (mm)');
set(gca,'FontSize',16 );

f4 = figure;
DR_1 = sqrt(DX_1.^2+DY_1.^2);
DR_2 = sqrt(DX_2.^2+DY_2.^2);
DR_3 = sqrt(DX_3.^2+DY_3.^2);
DR_rt_1 = sqrt(DX_rt_1.^2+DY_rt_1.^2);
DR_rt_2 = sqrt(DX_rt_2.^2+DY_rt_2.^2);
DR_rt_3 = sqrt(DX_rt_3.^2+DY_rt_3.^2);
##semilogy(100./(2*Rho(1,:)), sqrt(mean(DR_1.^2)));
##hold on;
##semilogy(100./(2*Rho(1,:)), sqrt(mean(DR_rt_1.^2)),'bo');
semilogy(100./(2*Rho(1,:)), sqrt(mean(DR_rt_1.^2))./sqrt(mean(DR_1.^2))-1,'b-');
hold on;
semilogy(100./(2*Rho(1,:)), sqrt(mean(DR_rt_2.^2))./sqrt(mean(DR_2.^2))-1,'r-');
semilogy(100./(2*Rho(1,:)), sqrt(mean(DR_rt_3.^2))./sqrt(mean(DR_3.^2))-1,'g-');

xlim([0 20]);ylim([0.001 1]);
xlabel('F number (-)');
ylabel('RMS error ratio (-)');
set(gca,'FontSize',16 );




##N_a = size(r_a2,2);
##N_b = size(r_b2,2);
##DR1 = zeros(N_a,N_b);
##DR2 = zeros(N_a,N_b);
##for i_a=1:1:N_a
##  disp([num2str(i_a) ' / ' num2str(N_a)]);
##  for i_b=1:1:N_b
##    lens1.r = [r_a1(i_a) r_a2(i_a) r_b1(i_b) r_b2(i_b)];                       
##    [aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, theta, EPD, conf);
##    DX_1 = aberdata.DX(5:end,:,1);
##    DY_1 = aberdata.DY(5:end,:,1);
##    DX_2 = aberdata.DX(5:end,:,end);
##    DY_2 = aberdata.DY(5:end,:,end);
##    
##    [N_X, N_Y] = size(DX_1(:,:,1));
##    DX_1_mu = mean(mean(DX_1));
##    DX_2_mu = mean(mean(DX_2));
##    DY_1_mu = mean(mean(DY_1));
##    DY_2_mu = mean(mean(DY_2));
##    
##    DX_1_sigma = sqrt(mean(mean((DX_1 - DX_1_mu).^2)));
##    DX_2_sigma = sqrt(mean(mean((DX_2 - DX_2_mu).^2)));
##    DY_1_sigma = sqrt(mean(mean((DY_1 - DY_1_mu).^2)));
##    DY_2_sigma = sqrt(mean(mean((DY_2 - DY_2_mu).^2)));
##    
##    DR1(i_a,i_b) = sqrt(DX_1_sigma^2 + DY_1_sigma^2);
##    DR2(i_a,i_b) = sqrt(DX_2_sigma^2 + DY_2_sigma^2);
##    
##  end
##end
##
##
##[val1, i_a_min1] = min(min(DR1,[],2),[],1);
##[val2, i_b_min1] = min(DR1(i_a_min1,:),[],2);
##[val3, i_a_min2] = min(min(DR2,[],2),[],1);
##[val4, i_b_min2] = min(DR2(i_a_min2,:),[],2);
##disp(['r_a1 ' 'r_a2 ' 'r_b1 ' 'r_b2']);
##disp([r_a1(i_a_min1) r_a2(i_a_min1) r_b1(i_b_min1) r_b2(i_b_min1)]);
##disp([r_a1(i_a_min2) r_a2(i_a_min2) r_b1(i_b_min2) r_b2(i_b_min2)]);
##val1
##val3
##
####r = [100 300];                % 曲率半径 (mm)
####b = [0 0];                    % 非球面係数
####d = [80 100];                 % 面間隔 (mm)
####n = [1.44 1];                 % 屈折率
####
####s_1 = -5000;                  % 第一面から物体面までの距離 (mm)
####t_1 = 50;               % 第一面から入射瞳までの距離 (mm)
####n_0 = 1;                      % 第一面から物体面における媒質の屈折率
####
####isFigLens = 'y';
####isFigAber = 'y';
####isRaytrace = 'y';
####
####
####%% ------------- データをバンドル ------------- %%
####lens1.r = r;                       
####lens1.b = b;           
####lens1.d = d;
####lens1.n = n;
####lens2.n_0 = n_0;
####lens2.s_1 = s_1;
####lens2.t_1 = t_1;
####conf.isFigLens = isFigLens;
####conf.isFigAber = isFigAber;
####conf.isRaytrace = isRaytrace;
####
####[aberdata, abercoef, gaussdata] = abercalc (lens1, lens2, phi, EPD, conf);