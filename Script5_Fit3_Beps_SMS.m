% Fit trivariate spline against the multiscale model data using B and
% epsilon as the variables. The script produces Figure 7 in the paper.

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
  load ./data/sms_data3
  
  % Interpolate H and lambda for uniform B
  Bref = Bx(:,1)';
  for i = 1 : length(sigxx)
    for j = 1 : length(sigxy)
      Hx2(:,i,j)    = interp1(Bx(:,i,j), Hx,           Bref, 'pchip', 'extrap');
      lamxx2(:,i,j) = interp1(Bx(:,i,j), lamxx(:,i,j), Bref, 'pchip', 'extrap');
      lamxy2(:,i,j) = interp1(Bx(:,i,j), lamxy(:,i,j), Bref, 'pchip', 'extrap');
    end;
  end;
  Bx = Bref;
  Hx = Hx2;
  lamxx = lamxx2;  
  lamxy = lamxy2; 
  clear Bref Hx2 lamxx2 lamxy2
  
  % Variables v' and w' which are not regularly distributed
  vp = (nu+1)/E*repmat(sigxx, length(Bx), 1, length(sigxx)) + 3/2*lamxx;
  wp = (nu+1)/E*repmat(reshape(sigxy,1,1,[]), length(Bx),length(sigxx),1) + lamxy;
  
  % Interpolate H, lambda and sigma for regular v and w
  v = linspace(min(vp(:)), max(vp(:)), length(sigxx));
  w = linspace(min(wp(:)), max(wp(:)), length(sigxy));
  for i = 1 : length(Bx)
    for j = 1 : length(sigxy)
      Hx2(i,:,j)    = interp1(vp(i,:,j), Hx(i,:,j),    v, 'pchip', 'extrap');
      sigxx2(i,:,j) = interp1(vp(i,:,j), sigxx,        v, 'pchip', 'extrap');
      sigxy2(i,:,j) = interp1(vp(i,:,j), sigxy,        v, 'pchip', 'extrap');
      lamxx2(i,:,j) = interp1(vp(i,:,j), lamxx(i,:,j), v, 'pchip', 'extrap');
      lamxy2(i,:,j) = interp1(vp(i,:,j), lamxy(i,:,j), v, 'pchip', 'extrap');
    end;
  end;
  for i = 1 : length(Bx)
    for j = 1 : length(sigxx)
      Hx2(i,j,:)    = interp1(squeeze(wp(i,j,:)), squeeze(Hx2(i,j,:)),    w, 'pchip', 'extrap');
      sigxx2(i,j,:) = interp1(squeeze(wp(i,j,:)), squeeze(sigxx2(i,j,:)), w, 'pchip', 'extrap');
      sigxy2(i,j,:) = interp1(squeeze(wp(i,j,:)), squeeze(sigxy2(i,j,:)), w, 'pchip', 'extrap');
      lamxx2(i,j,:) = interp1(squeeze(wp(i,j,:)), squeeze(lamxx2(i,j,:)), w, 'pchip', 'extrap');
      lamxy2(i,j,:) = interp1(squeeze(wp(i,j,:)), squeeze(lamxy2(i,j,:)), w, 'pchip', 'extrap');
    end;
  end;
  Hx = Hx2;
  sigxx = sigxx2;
  sigxy = sigxy2;
  lamxy = lamxy2;
  lamxx = lamxx2;
  clear vp wp Hx2 sigxx2 sigxy2 lamxx2 lamxy2
  
  % Auxiliary variables with scaling
  Bscale = max(Bx);
  exx_scale = max(abs(v));
  exy_scale = max(abs(w));
  u = Bx/Bscale;
  v = v/exx_scale;
  w = w/exy_scale;

  % Partial derivatives of the spline obtained from the multiscale model
  phi_u = Hx*Bscale;
  phi_v = -E/(1+nu)*lamxx*exx_scale;
  phi_w = -2*E/(1+nu)*lamxy*exy_scale;

%%% Fit spline
  
  tic
  s = fitSpline3(ordr, u, v, w, phi_u, phi_v, phi_w);
  toc
  
  % Keep scaling coefficients
  s.Bscale = Bscale;
  s.exx_scale = exx_scale;
  s.exy_scale = exy_scale;

  % Save results
  save ./splines/s3d_Beps s
  
%%% Plots
  
  % Calculate partial derivatives from the fitted spline
  sdu = fnval(fnder(s, [1 0 0]), {u,v,w});
  sdv = fnval(fnder(s, [0 1 0]), {u,v,w});
  sdw = fnval(fnder(s, [0 0 1]), {u,v,w});

  % Plot all data and show errors
  figure;
    hold on;
    plot(phi_u(:)/Bscale, 'b.-')
    plot(sdu(:)/Bscale, 'ro-')
    title(sprintf('H, error %g %%', 100*norm(phi_u(:)-sdu(:))/norm(phi_u(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(phi_v(:)/exx_scale, 'b.-')
    plot(sdv(:)/exx_scale, 'ro-')
    title(sprintf('\\sigma_{xx}, error %g %%', 100*norm(phi_v(:)-sdv(:))/norm(phi_v(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')
  figure;
    hold on;
    plot(phi_w(:)/exy_scale, 'b.-')
    plot(sdw(:)/exy_scale, 'ro-')
    title(sprintf('\\sigma_{xy}, error %g %%', 100*norm(phi_w(:)-sdw(:))/norm(phi_w(:))), 'FontSize', 14);
    legend('Multiscale', 'Spline')

   % 3-D plots with given indices of the B vector
   ibs = [3 6 15];
   [vv,ww] = ndgrid(v*exx_scale,w*exy_scale);
   vv = vv*1e6;
   ww = ww*1e6;
   for ii = 1 : length(ibs)
     ib = ibs(ii);

     figure(901);
       p2 = mesh(vv,ww,squeeze(sdu(ib,:,:))/Bscale);
       hold on;
       p1 = plot3(vv(:),ww(:),reshape(phi_u(ib,:,:),[],1)/Bscale, 'k.');
       xlabel('Variable {\itv} (ppm)', 'FontSize', 14);
       ylabel('Variable {\itw} (ppm)', 'FontSize', 14);
       zlabel('Field strength {\itH}_x (T)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(800,-1500,sdu(ib,end,1)/Bscale, ['{\itB} = ' num2str(Bx(ib), '%4.1f') ' T']);
     figure(902);
       p2 = mesh(vv,ww,-squeeze(sdv(ib,:,:))/exx_scale/1e6);
       hold on;
       p1 = plot3(vv(:),ww(:),-reshape(phi_v(ib,:,:),[],1)/exx_scale/1e6, 'k.');
       xlabel('Variable {\itv} (ppm)', 'FontSize', 14);
       ylabel('Variable {\itw} (ppm)', 'FontSize', 14);
       zlabel('Stress {\it\tau}_{xx} (MPa)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(800,-1500,-sdv(ib,end,1)/exx_scale/1e6, ['{\itB} = ' num2str(Bx(ib), '%4.1f') ' T']); 
     figure(903);
       p2 = mesh(vv,ww,-0.5*squeeze(sdw(ib,:,:))/exy_scale/1e6);
       hold on;
       p1 = plot3(vv(:),ww(:),-0.5*reshape(phi_w(ib,:,:),[],1)/exy_scale/1e6, 'k.');
       xlabel('Variable {\itv} (ppm)', 'FontSize', 14);
       ylabel('Variable {\itw} (ppm)', 'FontSize', 14);
       zlabel('Stress {\it\tau}_{xy} (MPa)', 'FontSize', 14);
       l = legend([p1 p2], 'Multiscale', 'Spline'); set(l, 'FontSize', 12);
       text(-800,-800,-0.5*sdw(ib,1,1)/exy_scale/1e6, ['{\itB} = ' num2str(Bx(ib), '%4.1f') ' T']); 
   end;
   
