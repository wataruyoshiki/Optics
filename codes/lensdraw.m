function retval = lensdraw(lens1, lens2, gaussdata, Yk, Zk, EPD, y_0, Y_obj, line)

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
  M = gaussdata.M;
  Md = gaussdata.Md;

  N = size(r,2);      % ���ܖʐ�
  lambda_0 = 1;
  Y_0 = (t_1-s_1)/(n_0*lambda_0)*y_0;                 % ���̍���


  %% ------------- �����Y�y�ь����̕`�� ------------- %%
  hold on;
  for i=1:1:N
    if(i==1)
    D = 0;              % ���_(��1)���瑪�肵����i�̋���(z����)
    else
      D = sum(d(1:i-1));
    end
    Y(i) = Yk(i);
    Z(i) = D + Zk(i);   % ���_(��1)���瑪�肵����i�̓��˓_�̋���(z����)

    %% ---- �����Y�`��`�� ---- %%
    y = (-1.1:0.01:1.1)*abs(Yk(i));
    z = sign(-r(i))*sqrt(r(i)^2-y.^2) + D + r(i);
    plot(z,y,'k-');    
  end
  
  %% ---- �����`�� ---- %%
  Y = [Y_0 Y Y_obj];
  Z = [s(1) Z sd(end)+sum(d(1:end-1))];
  if(abs(Z(1))>abs(Z(end)))
    Y(1) = (Y(2)-Y(1))/(Z(2)-Z(1))*(-Z(end)-Z(1)) + Y(1);
    Z(1) = -Z(end);
  end
  plot(Z,Y,line);

  %% ---- ����/�ˏo���`�� ---- %%
  Y_enp = [-EPD/2 EPD/2];
  Z_enp = [t(1) t(1)];
  Y_exp = [-Md*EPD/2 Md*EPD/2];
  Z_exp = [td(end)+sum(d(1:end-1)) td(end)+sum(d(1:end-1))];
  plot(Z_enp,Y_enp,'r-');
  plot(Z_exp,Y_exp,'g-');
  
  %% ---- �����`�� ---- %%
  plot([Z(1) Z(end)],[0 0],'k--');
  
  hold off;
  axis equal;

endfunction
