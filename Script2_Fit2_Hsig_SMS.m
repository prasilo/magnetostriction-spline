% Fit bivariate spline against the multiscale model data using H and sigma
% as the variables. The script produces Figure 4 in the paper.

clear all
close all
clc
addpath util

% Order of spline (4 means 3rd order)
ordr = 4;

%%% Load multiscale data

  % Data from multiscale model
  load ./data/sms_data2

  % Auxiliary variables with scaling
  Hscale = max(Hx);
  sscale = max(sigxx);
  u = Hx/Hscale;
  v = sigxx/sscale;

  % Partial derivatives of the spline obtained from the multiscale model
  phi_u = Bx*Hscale;
  phi_v = lamxx*sscale;

%%% Fit spline
  
  tic
  s = fitSpline2(ordr, u, v, phi_u, phi_v);
  toc
  
  % Keep scaling coefficients
  s.Hscale = Hscale;
  s.sscale = sscale;
  s.sx.Hscale = Hscale;
  s.sx.sscale = sscale;
  
  % Save results
  save ./splines/s2d_Hsig s
  
%%% Plots
  
  % Calculate partial derivatives from the fitted spline
  sdu = fnval(fnder(s, [1 0]), {u,v});
  sdv = fnval(fnder(s, [0 1]), {u,v});

  % Plot all data and show errors
  figure;
    hold on;
    plot(phi_u(:)/Hscale, 'b.-')
    plot(sdu(:)/Hscale, 'ro-')
    title(sprintf('B, error %g %%', 100*norm(phi_u(:)-sdu(:))/norm(phi_u(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(phi_v(:)/sscale, 'b.-')
    plot(sdv(:)/sscale, 'ro-')
    title(sprintf('{\\lambda}_{xx}, error %g %%', 100*norm(phi_v(:)-sdv(:))/norm(phi_v(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')

  % 3-D plots
  [hh,ss] = ndgrid(Hx, sigxx);
  figure;
    p2 = mesh(hh,ss/1e6,sdu/Hscale);
    hold on;
    p1 = plot3(hh(:),ss(:)/1e6,phi_u(:)/Hscale, 'k.');
    xlabel('Field strength {\itH} (A/m)', 'FontSize', 14);
    ylabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
    zlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
  figure;
    p2 = mesh(hh,ss/1e6,sdv/sscale*1e6);
    hold on;
    p1 = plot3(hh(:),ss(:)/1e6,phi_v(:)/sscale*1e6, 'k.');
    xlabel('Field strength {\itH} (A/m)', 'FontSize', 14);
    ylabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
    zlabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
    