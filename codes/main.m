clear all;
close all;

%% ------------- �p�����[�^��`  ------------- %%
r = [73.5017 -900.0000 9.0675 8.0000];                % �ȗ����a (mm)
b = [0 0 0 0];                    % �񋅖ʌW��
d = [0 96.939 0 54.286];                 % �ʊԊu (mm)
n = [1.5007 1 1.5007 1];                 % ���ܗ�

s_1 = -50000;                  % ���ʂ��畨�̖ʂ܂ł̋��� (mm)
t_1 = -0.0001;               % ���ʂ�����˓��܂ł̋��� (mm)
n_0 = 1;                      % ���ʂ��畨�̖ʂɂ�����}���̋��ܗ�
phi = deg2rad((0:1:10));      % ��p(���p) (rad.)
EPD = 50;                    % ���˓����a (mm)

isFigLens = 'y';
isFigAber = 'y';
isRaytrace = 'y';


%% ------------- �f�[�^���o���h�� ------------- %%
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