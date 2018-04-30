% Fit trivariate spline against the multiscale model data using H and sigma
% as the variables. The script produces Figure 5 in the paper.

clear all
close all
clc
addpath util

% Order of spline (4 means 3rd order)
ordr = 4;

%%% Load multiscale data

  % Data from multiscale model
  load ./data/sms_data3
  
  % Auxiliary variables with scaling
  Hscale = max(Hx);
  sxx_scale = max(sigxx);
  sxy_scale = max(sigxy);
  u = Hx/Hscale;
  v = sigxx/sxx_scale;
  w = sigxy/sxy_scale;

  % Partial derivatives of the spline obtained from the multiscale model
  phi_u = Bx*Hscale;
  phi_v = lamxx*sxx_scale;
  phi_w = 2*lamxy*sxy_scale;

%%% Fit spline
  
  tic
  s = fitSpline3(ordr, u, v, w, phi_u, phi_v, phi_w);
  toc

  % Keep scaling coefficients
  s.Hscale = Hscale;
  s.sxx_scale = sxx_scale;
  s.sxy_scale = sxy_scale;
  s.sx.Hscale = Hscale;
  s.sx.sxx_scale = sxx_scale;
  s.sx.sxy_scale = sxy_scale;
  
  % Save results
  save ./splines/s3d_Hsig s

%%% Plots  
  
  % Calculate partial derivatives from the fitted spline
  sdu = fnval(fnder(s, [1 0 0]), {u,v,w});
  sdv = fnval(fnder(s, [0 1 0]), {u,v,w});
  sdw = fnval(fnder(s, [0 0 1]), {u,v,w});
  
  % Plot all data and show errors
  figure;
    hold on;
    plot(phi_u(:)/Hscale, 'b.-')
    plot(sdu(:)/Hscale, 'ro-')
    title(sprintf('B, error %g %%', 100*norm(phi_u(:)-sdu(:))/norm(phi_u(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(phi_v(:)/sxx_scale, 'b.-')
    plot(sdv(:)/sxx_scale, 'ro-')
    title(sprintf('\\lambda_{xx}, error %g %%', 100*norm(phi_v(:)-sdv(:))/norm(phi_v(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(phi_w(:)/sxy_scale, 'b.-')
    plot(sdw(:)/sxy_scale, 'ro-')
    title(sprintf('\\lambda_{xy}, error %g %%', 100*norm(phi_w(:)-sdw(:))/norm(phi_w(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
    
   % 3-D plots with given indices of the H vector
   ihs = [3 8 21];
   [sxx,sxy] = ndgrid(sigxx, sigxy);
   for ii = 1 : length(ihs)
     ih = ihs(ii);

     figure(901);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdu(ih,:,:))/Hscale);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(phi_u(ih,:,:),[],1)/Hscale, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdu(ih,end,1)/Hscale, ['{\itH} = ' num2str(round(Hx(ih))) ' A/m']);
     figure(902);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdv(ih,:,:))/sxx_scale*1e6);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(phi_v(ih,:,:),[],1)/sxx_scale*1e6, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdv(ih,end,1)/sxx_scale*1e6, ['{\itH} = ' num2str(round(Hx(ih))) ' A/m']);
     figure(903);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdw(ih,:,:))/sxy_scale*1e6);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(phi_w(ih,:,:),[],1)/sxy_scale*1e6, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Magnetostriction {\it\lambda}_{xy} (ppm)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdw(ih,end,1)/sxy_scale*1e6, ['{\itH} = ' num2str(round(Hx(ih))) ' A/m']);
   end;
   