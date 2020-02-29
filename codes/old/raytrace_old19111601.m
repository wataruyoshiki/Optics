function [retval, Y, Z] = raytrace (r, b, d, n, s, sd, t, td, n_0, lambda_0, phi, EPD)

  eta_0 = EPD;
  y_0 = n_0*lambda_0*tan(phi);
  U_1 = -atan(lambda_0*eta_0/(t(1)-s(1)) - y_0/(n_0*lambda_0));
  L_1 = 1/tan(U_1)*(t(1)-s(1))/(n_0*lambda_0)*y_0 + s(1);

  N = size(r,2);
  L = zeros(1,N);
  Ld = zeros(1,N);
  U = zeros(1,N);
  Ud = zeros(1,N);
  I = zeros(1,N);
  Id = zeros(1,N);
  Y = zeros(1,N);
  Z = zeros(1,N);

  for i=1:1:N
    if(i==1)
      L(i) = L_1;
      U(i) = U_1;
      I(i) = asin((L(i)-r(i))/r(i)*sin(U(i)));
      Id(i) = asin(n_0/n(i)*sin(I(i)));
    else
      L(i) = Ld(i-1) - d(i-1);
      U(i) = Ud(i-1);
      I(i) = asin((L(i)-r(i))/r(i)*sin(U(i)));
      Id(i) = asin(n(i-1)/n(i)*sin(I(i)));
    end
    Ud(i) = U(i) + I(i) - Id(i);
    Ld(i) = sin(Id(i))/sin(Ud(i))*r(i) + r(i);
    Y(i) = r(i)*sin(U(i)+I(i));
    Z(i) = r(i) - r(i)*cos(U(i)+I(i));
  end
  
  retval = (Ld(end)-sd(end))*tan(Ud(end));

endfunction
