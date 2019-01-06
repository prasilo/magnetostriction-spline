% Fit bivariate spline against the multiscale model data using H and sigma
% as the variables. Some artificial numerical error is included in the
% results. The script produces Figure 8 (a) in the paper.

clear all
close all
clc
addpath util

% Numerical error for modifying the SMS model results
err = logspace(-5,0,21);

% Order of spline (4 means 3rd order)
ordr = 4;

%%% Load multiscale data

  % Loop over errors
  for ierr = 1 : length(err)

    % Data from multiscale model
    load ./data/sms_data2

    % Impose numerical into SMS model results
    Bx    = (1+err(ierr)/2)*Bx;
    lamxx = (1-err(ierr)/2)*lamxx;

    % Auxiliary variables with scaling
    Hscale = 1/max(abs(Bx(:)));
    sscale = 1/max(abs(lamxx(:)));
    u = Hx/Hscale;
    v = sigxx/sscale;

    % Partial derivatives of the spline obtained from the multiscale model
    psi_u = Bx*Hscale;
    psi_v = lamxx*sscale;

  %%% Fit spline

    tic
    s = fitSpline2(ordr, u, v, psi_u, psi_v);
    toc

  %%% Plots

    % Calculate partial derivatives from the fitted spline
    sdu = fnval(fnder(s, [1 0]), {u,v});
    sdv = fnval(fnder(s, [0 1]), {u,v});

    % Errors
    erru(ierr) = norm(psi_u(:)-sdu(:))/norm(psi_u(:));
    errv(ierr) = norm(psi_v(:)-sdv(:))/norm(psi_v(:));
    fprintf('Errors: \n');
    fprintf(' Bx:    %.3g %%\n', erru(ierr)*100);
    fprintf(' lamxx: %.3g %%\n', errv(ierr)*100);
  end
  
  % 3-D plots
  [hh,ss] = ndgrid(Hx, sigxx);
  figure;
    p2 = mesh(hh,ss/1e6,sdu/Hscale);
    hold on;
    p1 = plot3(hh(:),ss(:)/1e6,psi_u(:)/Hscale, 'k.');
    xlabel('Field strength {\itH}_x (A/m)', 'FontSize', 14);
    ylabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
    zlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
  figure;
    p2 = mesh(hh,ss/1e6,sdv/sscale*1e6);
    hold on;
    p1 = plot3(hh(:),ss(:)/1e6,psi_v(:)/sscale*1e6, 'k.');
    xlabel('Field strength {\itH}_x (A/m)', 'FontSize', 14);
    ylabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
    zlabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
    l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
  
  % Plot errors
  figure
  loglog(err*100, erru*100, 'x-'); hold on;
  loglog(err*100, errv*100, 'v-'); hold on;
  grid on
  xlabel('Input data error {\itr} (%)', 'FontSize', 14);
  ylabel('Fitting error {\itr}_{fit} (%)', 'FontSize', 14);
  set(gca, 'XTickLabel', num2str(10.^(-3:2)'))
  set(gca, 'YTickLabel', num2str(10.^(-3:2)'))
  l = legend('{\itB}_x', '{\it\lambda}_{xx}', 'Location', 'NorthWest'); set(l, 'FontSize', 12); 
    