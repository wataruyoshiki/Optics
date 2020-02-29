function gaussdata = gauss (lens1, lens2)

  %% ------------- �����Y�f�[�^���A���o���h�� ------------- %%
  r = lens1.r;                         
  b = lens1.b;           
  d = lens1.d;
  n = lens1.n;
  n_0 = lens2.n_0;
  s_1 = lens2.s_1;
  t_1 = lens2.t_1;

  N = size(r,2);                      % ���ܖʐ�

  %% ------------- �ϐ���`("���w�̌���"�Q��(d�̓_�b�V���̈�)) ------------- %%
  s = zeros(1,N);       % i�ʂ���̕��_����                     
  sd = zeros(1,N);      % i�ʂ���̑��_����
  t = zeros(1,N);       % i�ʂ���̓��˓�����
  td = zeros(1,N);      % i�ʂ���̎ˏo������
  h = zeros(1,N);       % i�ʂɂ����镨�̌����̒ʉߍ���
  k = zeros(1,N);       % ?�@�U�C�f�������v�Z�ɗp����
  K = zeros(1,N);       % i�ʂɂ����镨�̌�����Abbe�s�ϗ�
  L = zeros(1,N);       % i�ʂɂ����铵������Abbe�s�ϗ�

  
  %% ------------- �ߎ������ǐ� ------------- %%
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

  %% ------------- �o�̓f�[�^���o���h�� ------------- %%
  gaussdata.s = s;
  gaussdata.sd = sd;
  gaussdata.t = t;
  gaussdata.td = td;
  gaussdata.h = h;
  gaussdata.k = k;
  gaussdata.K = K;
  gaussdata.L = L;
  gaussdata.M = M;
  gaussdata.Md = Md;
  
  % [s, sd, t, td, h, k, K, L]
  
endfunction
