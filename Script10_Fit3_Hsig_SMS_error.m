% Fit trivariate spline against the multiscale model data using H and sigma
% as the variables.  Some artificial numerical error is included in the
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

  for ierr = 1 : length(err)

    % Data from multiscale model
    load ./data/sms_data3

    % Impose numerical into SMS model results
    Bx    = (1+err(ierr)/2)*Bx;
    lamxx = (1-err(ierr)/2)*lamxx;

    % Auxiliary variables with scaling
    Hscale = 1/max(abs(Bx(:)));
    sxx_scale = 1/max(abs(lamxx(:)));
    sxy_scale = 1/max(abs(lamxy(:)));
    u = Hx/Hscale;
    v = sigxx/sxx_scale;
    w = sigxy/sxy_scale;
    
    % Partial derivatives of the spline obtained from the multiscale model
    psi_u = Bx*Hscale;
    psi_v = lamxx*sxx_scale;
    psi_w = 2*lamxy*sxy_scale;

  %%% Fit spline

    tic
    s = fitSpline3(ordr, u, v, w, psi_u, psi_v, psi_w, 0);
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

    % Errors
    erru(ierr) = norm(psi_u(:)-sdu(:))/norm(psi_u(:));
    errv(ierr) = norm(psi_v(:)-sdv(:))/norm(psi_v(:));
    errw(ierr) = norm(psi_w(:)-sdw(:))/norm(psi_w(:));
    fprintf('Errors: \n');
    fprintf(' Bx:    %.3g %%\n', erru(ierr)*100);
    fprintf(' lamxx: %.3g %%\n', errv(ierr)*100);
    fprintf(' lamxy: %.3g %%\n', errw(ierr)*100);
  end
    
   % 3-D plots with given indices of the H vector
   ihs = [3 8 21];
   [sxx,sxy] = ndgrid(sigxx, sigxy);
   for ii = 1 : length(ihs)
     ih = ihs(ii);

     figure(901);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdu(ih,:,:))/Hscale);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(psi_u(ih,:,:),[],1)/Hscale, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Flux density {\itB}_x (T)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdu(ih,end,1)/Hscale, ['{\itH}_x = ' num2str(round(Hx(ih))) ' A/m']);
     figure(902);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdv(ih,:,:))/sxx_scale*1e6);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(psi_v(ih,:,:),[],1)/sxx_scale*1e6, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdv(ih,end,1)/sxx_scale*1e6, ['{\itH}_x = ' num2str(round(Hx(ih))) ' A/m']);
     figure(903);
       p2 = mesh(sxx/1e6,sxy/1e6,squeeze(sdw(ih,:,:))/sxy_scale*1e6);
       hold on;
       p1 = plot3(sxx(:)/1e6,sxy(:)/1e6,reshape(psi_w(ih,:,:),[],1)/sxy_scale*1e6, 'k.');
       xlabel('Stress {\it\sigma}_{xx} (MPa)', 'FontSize', 14);
       ylabel('Stress {\it\sigma}_{xy} (MPa)', 'FontSize', 14);
       zlabel('Magnetostriction {\it\lambda}_{xy} (ppm)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(105,-100,sdw(ih,end,1)/sxy_scale*1e6, ['{\itH}_x = ' num2str(round(Hx(ih))) ' A/m']);
   end
   
  % Plot errors
  figure
  loglog(err*100, erru*100, 'x-'); hold on;
  loglog(err*100, errv*100, 'v-'); hold on;
  loglog(err*100, errw*100, 'o-'); hold on;
  grid on
  xlabel('Input data error {\itr} (%)', 'FontSize', 14);
  ylabel('Fitting error {\itr}_{fit} (%)', 'FontSize', 14);
  set(gca, 'XTickLabel', num2str(10.^(-3:2)'))
  set(gca, 'YTickLabel', num2str(10.^(-3:2)'))
  l = legend('{\itB}_x', '{\it\lambda}_{xx}', '{\it\lambda}_{xy}', 'Location', 'NorthWest'); set(l, 'FontSize', 12); 
   
   