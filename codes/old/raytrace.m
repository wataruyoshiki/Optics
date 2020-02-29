function [retval, Y, Z] = raytrace (lens1, lens2, gaussdata, y_0, rho)

  %% ------------- �����Y�f�[�^���A���o���h�� ------------- %%
  r = lens1.r;                         
  b = lens1.b;           
  d = lens1.d;
  n = lens1.n;
  n_0 = lens2.n_0;
  s_1 = lens2.s_1;
  t_1 = lens2.t_1;
  
  %% ------------- �K�E�X���w�v�Z�f�[�^���A���o���h�� ------------- %%
  s = gaussdata.s;
  sd = gaussdata.sd;
  t = gaussdata.t;
  td = gaussdata.td;

  lambda_0 = 1;
  eta_0 = rho';
  U_1 = -atan(lambda_0.*eta_0/(t(1)-s(1)) - y_0/(n_0*lambda_0));
  L_1 = 1./tan(U_1)*(t(1)-s(1))/(n_0*lambda_0)*y_0 + s(1);

  N2 = size(r,2);           % �����Y�ʐ�
  N1 = size(rho',1);         % ���a�̏�����
  L = zeros(N1,N2);         % ��i�ւ̓��ˌ����������ƌ���鋗��
  Ld = zeros(N1,N2);        % ��i����̋��܌����������ƌ���鋗��
  U = zeros(N1,N2);         % ��i�ւ̓��ˌ����������ƂȂ��p
  Ud = zeros(N1,N2);        % ��i����̋��܌����������ƂȂ��p
  I = zeros(N1,N2);         % ��i�ւ̓��ˊp
  Id = zeros(N1,N2);        % ��i����̋��܊p
  Y = zeros(N1,N2);         % ��i�ւ̓��ˍ���
  Z = zeros(N1,N2);         % ��i�̓��ˍ��W�̃U�O��

  %% ------------- �������ǐ� ------------- %%
  for i=1:1:N2
    if(i==1)
      L(:,i) = L_1;
      U(:,i) = U_1;
      I(:,i) = asin((L(:,i)-r(i))./r(i).*sin(U(:,i)));
      Id(:,i) = asin(n_0/n(i).*sin(I(:,i)));
    else
      L(:,i) = Ld(:,i-1) - d(i-1);
      U(:,i) = Ud(:,i-1);
      I(:,i) = asin((L(:,i)-r(i))/r(i).*sin(U(:,i)));
      Id(:,i) = asin(n(i-1)/n(i).*sin(I(:,i)));
    end
    Ud(:,i) = U(:,i) + I(:,i) - Id(:,i);
    Ld(:,i) = sin(Id(:,i))./sin(Ud(:,i))*r(i) + r(i);
    Y(:,i) = r(i).*sin(U(:,i)+I(:,i));
    Z(:,i) = r(i) - r(i)*cos(U(:,i)+I(:,i));
  end
  
  %% ------------- ���ʂɂ������������ ------------- %%
  retval = (Ld(:,end)-sd(:,end)).*tan(Ud(:,end));

endfunction
