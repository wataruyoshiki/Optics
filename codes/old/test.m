% �R�[�h�̈ꕔ���֐���
% 3���������ǐՎ���

clear all;
close all;

%% ------------- �p�����[�^��`  ------------- %%
r = [100 500];                % �ȗ����a (mm)
b = [0 0];                    % �񋅖ʌW��
d = [50 100];                 % �ʊԊu (mm)
n = [1.44 1];                 % ���ܗ�

s_1 = -500;                   % ���ʂ��畨�̖ʂ܂ł̋��� (mm)
t_1 = -0.00001;               % ���ʂ�����˓��܂ł̋��� (mm)
n_0 = 1;                      % ���ʂ��畨�̖ʂɂ�����}���̋��ܗ�
phi = deg2rad(5);             % ��p(���p) (rad.)
EPD = 100;                    % ���˓����a (mm)


%% ------------- �ϐ���` ------------- %%
lambda_0 = 1;                       % 1�ȊO�ɐݒ肵�Ȃ����ƁB
y_0 = n_0*lambda_0*tan(phi);        % ���K�����ꂽ���̖ʂɂ����镨�_����
Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0; % ���̍���
N = size(r,2);                      % ���ܖʐ�

s = zeros(1,N);               % ���w�̌����Q�� (d�̓_�b�V���̈�)
sd = zeros(1,N);
t = zeros(1,N);
td = zeros(1,N);
h = zeros(1,N);
k = zeros(1,N);
K = zeros(1,N);
L = zeros(1,N);

B = zeros(1,N);               % ���ʎ����W��
C = zeros(1,N);               % ��_�����W��
D = zeros(1,N);               % ���ʘp�Ȏ����W��
E = zeros(1,N);               % �c�Ȏ����W��
F = zeros(1,N);               % �R�}�����W��


%% ------------- �K�E�X���w�v�Z ------------- %%
for i=1:1:N
  if(i==1)
    s(i) = s_1;
    sd(i) = ((1-n_0/n(i))*1/r(i) + n_0/n(i)*1/s(i))^(-1);
    t(i) = t_1;
    td(i) = ((1-n_0/n(i))*1/r(i) + n_0/n(i)*1/t(i))^(-1);
    h(i) = s_1/(t_1-s_1);
    k(i) = t_1*(t_1-s_1)/(n_0*s_1);
  else
    s(i) = sd(i-1)-d(i-1);
    sd(i) = ((1-n(i-1)/n(i))*1/r(i) + n(i-1)/n(i)*1/s(i))^(-1);
    t(i) = td(i-1)-d(i-1);
    td(i) = ((1-n(i-1)/n(i))*1/r(i) + n(i-1)/n(i)*1/t(i))^(-1);
    h(i) = h(i-1)*s(i)/sd(i-1);
    k_tmp = 0;
    for j=1:1:i-1
      k_tmp = k_tmp + d(j)/(n(j)*h(j)*h(j+1));
    end
    k(i) = k(1) + k_tmp;
  end
  K(i) = n(i)*(1/r(i) - 1/sd(i));
  L(i) = n(i)*(1/r(i) - 1/td(i));
end

M = n_0*prod(n(1:end-1))*prod(sd)/(prod(n)*prod(s));  % ���{��
Md = n_0*prod(n(1:end-1))*prod(td)/(prod(n)*prod(t)); % ���{��


%% ------------- ���H�}�`�� ------------- %%
[tmp, Yk, Zk] = raytrace(r, b, d, n, s, sd, t, td, n_0, lambda_0, phi, EPD/2);
lensdraw(r, b, d, n, s, sd, t, td, Yk, Zk, EPD, Md);


%% ------------- �����W���v�Z ------------- %%
for i=1:1:N
  if(i==1)
    N1 = n(1) - n_0;
    N2 = 1/n(1)^2 - 1/n_0^2;
    NS = 1/(n(1)*sd(1)) - 1/(n_0*s(1));
  else
    N1 = n(i) - n(i-1);
    N2 = 1/n(i)^2 - 1/n(i-1)^2;
    NS = 1/(n(i)*sd(i)) - 1/(n(i-1)*s(i));
  end
  B(i) = 1/2*( h(i)^4*(b(i)/r(i)^3)*N1 + h(i)^4*K(i)^2*NS );
  C(i) = 1/2*( h(i)^4*k(i)^2*(b(i)/r(i)^3)*N1 + (1 + h(i)^2*k(i)*K(i))^2*NS );
  D(i) = 1/2*( h(i)^4*k(i)^2*(b(i)/r(i)^3)*N1 + h(i)^2*k(i)*K(i)*(2 + h(i)^2*k(i)*K(i))*NS - K(i)*N2 );
  E(i) = 1/2*( h(i)^4*k(i)^3*(b(i)/r(i)^3)*N1 + k(i)*(1 + h(i)^2*k(i)*K(i))*(2 + h(i)^2*k(i)*K(i))*NS - (1 + h(i)^2*k(i)*K(i))/h(i)^2*N2 );
  F(i) = 1/2*( h(i)^4*k(i)*(b(i)/r(i)^3)*N1 + h(i)^2*K(i)*(1 + h(i)^2*k(i)*K(i))*NS );
end
B_tot = sum(B);     % �����W���̊e�ʍ��v
C_tot = sum(C);
D_tot = sum(D);
E_tot = sum(E);
F_tot = sum(F);


%% ------------- ���Z ------------- %%
P = (n(1)-1)*(1/r(1) - 1/r(2));
K = -1/s(1) - P/2;
sigma = (n(1)-1)*(1/r(1) + 1/r(2));
beta = 0;

U = 1/2*beta + n(1)^2/(8*(n(1)-1)^2)*P^3 - n(1)/(2*(n(1)+2))*K^2*P + 1/(2*n(1)*(n(1)+2))*P*((n(1)+2)/(2*(n(1)-1))*sigma + 2*(n(1)+1)*K)^2;
V = 1/(2*n(1))*P*((n(1)+1)/(2*(n(1)-1))*sigma + (2*n(1)+1)*K);
B_test = h(1)^4*U;
F_test = h(1)^4*k(1)*U + h(1)^2*V;
C_test = h(1)^4*k(1)^2*U + 2*h(1)^2*k(1)*V + 1/2*P;
D_test = h(1)^4*k(1)^2*U + 2*h(1)^2*k(1)*V + (n(1)+1)/(2*n(1))*P;
E_test = h(1)^4*k(1)^3*U + 3*h(1)^2*k(1)^2*V + k(1)*(3*n(1)+1)/(2*n(1))*P;

data_test = [
  B_tot B_test;
  C_tot C_test;
  D_tot D_test;
  E_tot E_test;
  F_tot F_test
  ];
disp(data_test);


%% ------------- �������ʌv�Z ------------- %%
rho = ((0:0.01:1)+0.000001)*(EPD/2);       % �ˏo���ɂ�����a�������� (mm)
theta = deg2rad(0:5:360)';
[Rho, Theta] = meshgrid(rho,theta);
Dx = B_tot*Rho.^2.*sin(Theta) - 2*F_tot*y_0*Rho.^2.*sin(Theta).*cos(Theta) + D_tot*y_0^2*Rho.*sin(Theta);
Dy = B_tot*Rho.^3.*cos(Theta) - F_tot*y_0*Rho.^2.*(1+2*cos(Theta).^2) + (2*C_tot+D_tot)*y_0^2*Rho.*cos(Theta) - E_tot*y_0^3;
DX = (td(end)-sd(end))/(n(end)*Md)*Dx;
DY = (td(end)-sd(end))/(n(end)*Md)*Dy;


%% ------------- �������ǐ� ------------- %%
N_rho = size(rho,2);
Y_rt = zeros(1,N_rho);
for i=1:1:N_rho
  Y_rt(i) = raytrace (r, b, d, n, s, sd, t, td, n_0, lambda_0, phi, rho(i));
end

figure;plot(rho,Y_rt-Y_0*M,'r-',rho,DY(1,:),'b-');