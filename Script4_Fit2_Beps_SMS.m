% Fit bivariate spline against the multiscale model data using B and
% epsilon as the variables. The script produces Figure 6 in the paper.

clear all
close all
clc
addpath util

% Young modulus and Poisson ratio
nu = 0.34;
E = 183e9;

% Order of spline (4 means 3rd order)
ordr = 4;

%%% Load multiscale data

  % Data from multiscale model
  load ./data/sms_data2

  % Interpolate H and lambda for uniform B
  Bref = Bx(:,1)';
  for i = 1 : length(sigxx)
    Hx2(:,i)    = interp1(Bx(:,i), Hx,         Bref, 'pchip', 'extrap');
    lamxx2(:,i) = interp1(Bx(:,i), lamxx(:,i), Bref, 'pchip', 'extrap');
  end
  Bx = Bref;
  Hx = Hx2;
  lamxx = lamxx2;
  clear Bref Hx2 lamxx2
  
  % Variable v' which is not regularly distributed
  vp = (nu+1)/E*repmat(sigxx, length(Bx), 1) + 3/2*lamxx;
  
  % Interpolate H, lambda and sigma for regular v
  v = linspace(min(vp(:)), max(vp(:)), length(sigxx));
  for i = 1 : length(Bx)
    Hx2(i,:)    = interp1(vp(i,:), Hx(i,:),    v, 'pchip', 'extrap');
    sigxx2(i,:) = interp1(vp(i,:), sigxx,      v, 'pchip', 'extrap');
    lamxx2(i,:) = interp1(vp(i,:), lamxx(i,:), v, 'pchip', 'extrap');
  end
  Hx = Hx2;
  sigxx = sigxx2;
  lamxx = lamxx2;
  clear vp sigxx2 lamxx2 Hx2

  % Auxiliary variables with scaling
  Bscale = 1/max(abs(Hx(:)));
  escale = 1/max(lamxx(:)*E/(1+nu));
  u = Bx/Bscale;
  v = v/escale;

  % Partial derivatives of the spline obtained from the multiscale model
  psi_u = Hx*Bscale;
  psi_v = -E/(1+nu)*lamxx*escale;
  
%%% Fit spline

  tic
  s = fitSpline2(ordr, u, v, psi_u, psi_v, 1);
  toc
  
  % Keep scaling coefficients
  s.Bscale = Bscale;
  s.escale = escale;
  s.sx.Bscale = Bscale;
  s.sx.escale = escale;
  save ./splines/s2d_Beps s
  
%%% Plots  
  
  % Calculate partial derivatives from the fitted spline
  sdu = fnval(fnder(s, [1 0]), {u,v});
  sdv = fnval(fnder(s, [0 1]), {u,v});

  % Errors
  erru = norm(psi_u(:)-sdu(:))/norm(psi_u(:));
  errv = norm(psi_v(:)-sdv(:))/norm(psi_v(:));
  fprintf('Errors: \n');
  fprintf(' Hx:    %.3g %%\n', erru*100);
  fprintf(' tauxx: %.3g %%\n', errv*100);

  % Plot all data and show errors
  figure;
    hold on;
    plot(psi_u(:)/Bscale, 'b.-')
    plot(sdu(:)/Bscale, 'ro-')
    title(sprintf('H, error %g %%', 100*erru), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(psi_v(:)/escale, 'b.-')
    plot(sdv(:)/escale, 'ro-')
    title(sprintf('\\sigma_{xx}, error %g %%', 100*errv), 'FontSize', 14);
    legend('Multiscale', 'Spline')
    
  % 3-D plots
  [bb,ee] = ndgrid(Bx, v*escale);
  figure;
    p2 = mesh(bb,ee*1e6,sdu/Bscale);
    hold on;
    p1 = plot3(bb(:),ee(:)*1e6,psi_u(:)/Bscale, 'k.');
    xlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
    ylabel('Variable {\itv} (ppm)', 'FontSize', 14);
    zlabel('Field strength {\itH}_x (A/m)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
  figure;
    p2 = mesh(bb,ee*1e6, -sdv/escale/1e6);
    hold on;
    p1 = plot3(bb(:),ee(:)*1e6, -psi_v(:)/escale/1e6, 'k.');
    xlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
    ylabel('Variable {\itv} (ppm)', 'FontSize', 14);
    zlabel('Stress {\it\tau}_{xx} (MPa)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
    