% Fit bivariate spline against measurement data using B and sigma as the
% variables. The script produces Figure 10 in the paper.

clear all
close all
clc
addpath util

% Young modulus and Poisson ratio
nu = 0.34;
E = 183e9;

% Order of spline (4 means 3rd order)
ordr = 4;

% Which stresses values to plot. The available values are [-35 -30 -20 -15
% -10 -5  0  5  10  15  40  50  60  70  80  90 100] MPa
sref = [-35 0 15 80];

%%% Load measurement data

  % Measurement data
  load ./data/meas_data2.mat
  ind = find(ismember(sigxx/1e6, sref));
  
  % Plot measurements
  num2strcell = @(a, pre, ext) reshape(cellstr([repmat(pre, [numel(a) 1])  num2str(a(:)) repmat(ext, [numel(a) 1])]), size(a));
  style = {'x-', 'v-', 'o-', 's-'};
  for i = 1 : length(ind)
    figure(1);
      plot(Hx(1:5:end,ind(i))', Bx(1:5:end)', style{i}); hold on;
      xlabel('Field strength {\itH} (A/m)', 'FontSize', 14);
      ylabel('Flux density {\itB}_x (T)', 'FontSize', 14);
    figure(2);
      plot(Hx(1:5:end,ind(i))', lamxx(1:5:end,ind(i))'*1e6, style{i}); hold on;
      xlabel('Field strength {\itH} (A/m)', 'FontSize', 14);
      ylabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
  end;
  figure(1); leg = legend(num2strcell(sref, '{\it\sigma}_{xx} =', ' MPa'), 'Location', 'SouthEast'); set(leg, 'FontSize', 12); legend boxoff
  drawnow;

  % Variable v' which is not regularly distributed
  vp = (nu+1)/E*repmat(sigxx, length(Bx), 1) + 3/2*lamxx;

  % Keep origindal measurement data
  Hx0 = Hx;
  sigxx0 = sigxx;
  lamxx0 = lamxx;
  
  % Interpolate H, lambda and sigma for regular v
  v = linspace(min(vp(:)), max(vp(:)), length(sigxx));
  for i = 1 : length(Bx)
    Hx2(i,:)    = interp1(vp(i,:), Hx(i,:),    v, 'pchip', 'extrap');
    sigxx2(i,:) = interp1(vp(i,:), sigxx,      v, 'pchip', 'extrap');
    lamxx2(i,:) = interp1(vp(i,:), lamxx(i,:), v, 'pchip', 'extrap');
  end;
  Hx = Hx2;
  sigxx = sigxx2;
  lamxx = lamxx2;
  clear vp sigxx2 lamxx2 Hx2

  % Auxiliary variables with scaling
  Bscale = max(Bx);
  escale = max(abs(v));
  u = Bx/Bscale;
  v = v/escale;

  % Stress caused by magnetostriction
  tauxx = -E/(1+nu)*lamxx;
  
  % Partial derivatives of the spline obtained from the multiscale model
  psi_u = Hx*Bscale;
  psi_v = tauxx*escale;

%%% Fit spline

  tic
  s = fitSpline2(4, u, v, psi_u, psi_v);
  toc
  
%%% Plots
  
  % Calculate partial derivatives from the fitted spline
  sdu = fnval(fnder(s, [1 0]), {u,v});
  sdv = fnval(fnder(s, [0 1]), {u,v});

  % Plot all data and show errors
  figure;
    hold on;
    plot(psi_u(:)/Bscale, 'b.-')
    plot(sdu(:)/Bscale, 'ro-')
    title(sprintf('B, error %g %%', 100*norm(psi_u(:)-sdu(:))/norm(psi_u(:))), 'FontSize', 14);
    legend('Measured', 'Spline')
  figure;
    hold on;
    plot(psi_v(:)/escale, 'b.-')
    plot(sdv(:)/escale, 'ro-')
    title(sprintf('{\\tau}_{xx}, error %g %%', 100*norm(psi_v(:)-sdv(:))/norm(psi_v(:))), 'FontSize', 14);
    legend('Measured', 'Spline')

  % Plot B(H) and lambda(B) curves with more points in B to see possible
  % oscillations in the splines
  style = {'x', 'v', 'o', 's'};
  for i = 1 : length(ind)
    
    % Auxiliary variables
    uu = Bx/Bscale;
    vv = ((nu+1)/E*repmat(sigxx0(ind(i)), length(Bx), 1) + 3/2*lamxx0(:,ind(i)))/escale;

    % Partial derivatives from spline
    sdu = fnval(fnder(s, [1 0]), [uu(:)'; vv(:)']);
    sdv = fnval(fnder(s, [0 1]), [uu(:)'; vv(:)']);
    
    figure(11);
    set(gca, 'ColorOrderIndex', i);
      p1(i) = plot(Hx0(1:5:end,ind(i)), Bx(1:5:end), style{i}); hold on;
      plot(sdu/Bscale, uu*Bscale, '-', 'Color', get(p1(i), 'color')); hold on;
      xlabel('Field strength {\itH}_x (A/m)', 'FontSize', 14);
      ylabel('Flux density {\itB} (T)', 'FontSize', 14);
    figure(12)
    set(gca, 'ColorOrderIndex', i);
      p12 = plot(Bx(1:5:end)', lamxx0(1:5:end,ind(i))*1e6, style{i}); hold on;
      plot(uu*Bscale, -(1+nu)/E*sdv/escale*1e6, '-', 'Color', get(p12, 'color')); hold on;
      xlabel('Flux density {\itB} (T)', 'FontSize', 14);
      ylabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
  end;    

  % Set legend
  figure(11)
  p2 = plot(0,0,'k-');
  leg = legend([p1 p2], [num2strcell(sref, 'Measured {\it\sigma}_{xx} =', ' MPa'), 'Spline'], 'Location', 'SouthEast'); set(leg, 'FontSize', 12); legend boxoff
  xlim([0 3000])
