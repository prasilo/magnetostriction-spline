% Compare spline-based invariant model and simplified multiscale model
% under shear stress. The script produces Figure 10 in the paper.

close all;
clear all;
clc;
addpath util util/functions_Hsig

mu0 = 4*pi*1e-7;

% Indices for changing from Voigt to tensor notation and vice versa
ivt = [1 6 5; 6 2 4; 5 4 3];
ivv = [1 5 9 8 7 4];

% H vector
H = [1077; 0; 0];

% Loop through bi-variate and trivariate cases
for flag3 = 0 : 1

  % Load splines
  if flag3 == 0
    load ./splines/s2d_Hsig;
    style = '--';
  else
    load ./splines/s3d_Hsig;
    style = '-';
  end

  % Starting stress
  sig0 = [0 1 0; 1 0 0; 0 0 0]*50e6;

  % Angles
  theta = linspace(0,pi,21);

  % Rotation matrix
  T = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];

  % Run models
  % Loop field strengths
  for i = 1 : length(theta)
    ss = T(theta(i))*sig0*T(theta(i))';

    % Evaluate multiscale model
    [bb,ll] = smsfunc(H, ss);

    % Relative permeability
    Bx_ms(i) = bb(1);
    By_ms(i) = bb(2);
    lxx_ms(i) = ll(1,1);
    lyy_ms(i) = ll(2,2);
    lzz_ms(i) = ll(3,3);
    lxy_ms(i) = ll(1,2);
    lyz_ms(i) = ll(2,3);
    lzx_ms(i) = ll(3,1);

    % Run invariant model
    [bb,ll]= ifunc_Hsig(s, H, ss(ivv)');
    Bx_i(i) = bb(1);
    By_i(i) = bb(2);
    lxx_i(i) = ll(1);
    lyy_i(i) = ll(2);
    lzz_i(i) = ll(3);
    lxy_i(i) = ll(6);
    lyz_i(i) = ll(4);
    lzx_i(i) = ll(5);

    % Run invariant model
    [bb,ll]= ifunc_Hsig(s.sx, H, ss(ivv)');
    Bx_i_xt(i) = bb(1);
    By_i_xt(i) = bb(2);
    lxx_i_xt(i) = ll(1);
    lyy_i_xt(i) = ll(2);
    lzz_i_xt(i) = ll(3);
    lxy_i_xt(i) = ll(6);
    lyz_i_xt(i) = ll(4);
    lzx_i_xt(i) = ll(5);
  end

  % Plots
  figure(1);
    hold on;
    set(gca, 'ColorOrderIndex', 1);
    px1 = plot(theta*180/pi,Bx_ms, 'x');
    plot(theta/pi*180,Bx_i, style, 'Color', get(px1, 'Color'));
    ps1(flag3+1) = plot(0,0,['k' style]);
    set(gca, 'ColorOrderIndex', 2);
    py1 = plot(theta*180/pi,By_ms, 'v');
    plot(theta/pi*180,By_i, style, 'Color', get(py1, 'Color'));
    xlabel('Angle {\it\theta} (°)', 'FontSize', 14);  
    ylabel('Flux density (T)', 'FontSize', 14);  
  figure(2);
    hold on;
    set(gca, 'ColorOrderIndex', 1);
    px2 = plot(theta*180/pi,lxx_ms*1e6, 'x');
    plot(theta/pi*180,lxx_i*1e6, style, 'Color', get(px2, 'Color'));
    ps2(flag3+1) = plot(0,0,['k' style]);
    set(gca, 'ColorOrderIndex', 2);
    py2 = plot(theta*180/pi,lyy_ms*1e6, 'v');
    plot(theta/pi*180,lyy_i*1e6, style, 'Color', get(py2, 'Color'));
    xlabel('Angle {\it\theta} (°)', 'FontSize', 14);  
    ylabel('Magnetostriction (ppm)', 'FontSize', 14);  
  figure(3);
    hold on;
    set(gca, 'ColorOrderIndex', 1);
    px3 = plot(theta*180/pi,lxy_ms*1e6, 'x');
    plot(theta/pi*180,lxy_i*1e6, style, 'Color', get(px3, 'Color'));
    ps3(flag3+1) = plot(0,0,['k' style]);
    set(gca, 'ColorOrderIndex', 2);
    py3 = plot(theta*180/pi,lzz_ms*1e6, 'v');
    plot(theta/pi*180,lzz_i*1e6, style, 'Color', get(py3, 'Color'));
    xlabel('Angle {\it\theta} (°)', 'FontSize', 14);  
    ylabel('Magnetostriction (ppm)', 'FontSize', 14);  
  drawnow;
end

% Set legends etc.
figure(1)
  l = legend([px1 py1 ps1], '{\itB}_x multiscale', '{\itB}_y multiscale','Bivariate spline', 'Trivariate spline'); set(l, 'FontSize', 12)
figure(2)
  l = legend([px2 py2 ps2], '{\it\lambda}_{xx} multiscale', '{\it\lambda}_{yy} multiscale','Bivariate spline', 'Trivariate spline'); set(l, 'FontSize', 12)
figure(3)
  l = legend([px3 py3 ps3], '{\it\lambda}_{xy} multiscale', '{\it\lambda}_{zz} multiscale','Bivariate spline', 'Trivariate spline'); set(l, 'FontSize', 12)
