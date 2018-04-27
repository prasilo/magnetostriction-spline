% Compare spline-based invariant model and simplified multiscale model
% under different multiaxial stresses. The script produces Figure 8 in the
% paper.

clear all
close all
clc
addpath util util/functions_Hsig

mu0 = 4*pi*1e-7;             

% Indices for changing from Voigt to tensor notation
ivt = [1 6 5; 6 2 4; 5 4 3];

% H vector and stress values to loop through
H = [1077; 0; 0];
sval = linspace(-100,100,51)*1e6;

% Loop through bi-variate and trivariate cases
for flag3 = 0 : 1

  % Load splines
  if flag3 == 0
    close all;
    load ./splines/s2d_Hsig;
  else
    load ./splines/s3d_Hsig;
  end;

  % Different multiaxial stresses
  sig = {}; slab = {}; mar = {};
  sig{end+1} = [1 0 0 0 0 0];  slab{end+1} = 'Uniaxial';    mar{end+1} = 'x';
  sig{end+1} = [1 1 0 0 0 0];  slab{end+1} = 'Equibiaxial'; mar{end+1} = 'v';
  sig{end+1} = [1 1 1 0 0 0];  slab{end+1} = 'Hydrostatic'; mar{end+1} = 'o';
  sig{end+1} = [1 -1 0 0 0 0]; slab{end+1} = 'Pure shear';  mar{end+1} = 's';

  % Loop thourgh each case
  p2 = [];
  slab2 = {};
  for isig = 1 : length(sig)

    % Run multiscale model (this takes some time)
    fprintf('Calculating multiscale model results...\n')
    for is = 1:length(sval)
      fprintf('%d / %d\n', is, length(sval));

      % Stress tensor
      ss = sig{isig}(ivt)*sval(is);

      % Evaluate multiscale model
      [bb,ll] = smsfunc(H, ss);

      % Relative permeability
      mur_ms(is,isig) = bb(1)/H(1)/mu0;
      lxx_ms(is,isig) = ll(1,1);
      lyy_ms(is,isig) = ll(2,2);
      lzz_ms(is,isig) = ll(3,3);
    end;

    % Run invariant model
    [bb,ll]= ifunc_Hsig(s, H*ones(size(sval)), sig{isig}'*sval);
    mur_i(:,isig) = bb(1,:)./H(1)/mu0;
    lxx_i(:,isig) = ll(1,:);
    lyy_i(:,isig) = ll(2,:);
    lzz_i(:,isig) = ll(3,:);

    % Run invariant model with extrapolation
    [bb,ll]= ifunc_Hsig(s.sx, H*ones(size(sval)), sig{isig}'*sval);
    mur_i_xt(:,isig) = bb(1,:)./H(1)/mu0;
    lxx_i_xt(:,isig) = ll(1,:);
    lyy_i_xt(:,isig) = ll(2,:);
    lzz_i_xt(:,isig) = ll(3,:);

    % Volume changes in magnetostriction from multiscale model and
    % invariant model
    V_ms = norm(lxx_ms + lyy_ms + lzz_ms);
    V_i = norm(lxx_i + lyy_i + lzz_i);
    V_i_xt = norm(lxx_i_xt + lyy_i_xt + lzz_i_xt);

    % Permeability plot
    figure(1);
      set(gca, 'ColorOrderIndex', isig);
      hold on;
      p(isig) = plot(sval(1:5:end)/1e6,mur_ms(1:5:end,isig), mar{isig});
      plot(sval/1e6,mur_i(:,isig),'-', 'Color', get(p(isig), 'Color'))
      plot(sval/1e6,mur_i_xt(:,isig),'--', 'Color', get(p(isig), 'Color'))

    % Magnetostriction plot
    figure(2);
      set(gca, 'ColorOrderIndex', isig);
      hold on;
      plot(sval(1:5:end)/1e6,lxx_ms(1:5:end,isig)*1e6, mar{isig});
      plot(sval/1e6,lxx_i(:,isig)*1e6,'-', 'Color', get(p(isig), 'Color'))
      plot(sval/1e6,lxx_i_xt(:,isig)*1e6,'--', 'Color', get(p(isig), 'Color'))
    drawnow;
     
    % Magnetostriction yy and zz plots for equibiaxial and pure shear cases
    if strcmp(slab{isig}, 'Equibiaxial') || strcmp(slab{isig}, 'Pure shear')
      figure(3);
        set(gca, 'ColorOrderIndex', isig);
        hold on;
        p2(end+1) = plot(sval(1:5:end)/1e6,lyy_ms(1:5:end,isig)*1e6, mar{isig});
        slab2{end+1} = slab{isig};
        plot(sval/1e6,lyy_i(:,isig)*1e6,'-', 'Color', get(p(isig), 'Color'))
        plot(sval/1e6,lyy_i_xt(:,isig)*1e6,'--', 'Color', get(p(isig), 'Color'))
      figure(4);
        set(gca, 'ColorOrderIndex', isig);
        hold on;
        plot(sval(1:5:end)/1e6,lzz_ms(1:5:end,isig)*1e6, mar{isig});
        plot(sval/1e6,lzz_i(:,isig)*1e6,'-', 'Color', get(p(isig), 'Color'))
        plot(sval/1e6,lzz_i_xt(:,isig)*1e6,'--', 'Color', get(p(isig), 'Color'))
      drawnow;
    end;
  end;
end;

% Set legends etc.
figure(1);
  p(end+1) = plot(sval(end)/1e6,mur_ms(end,end,end),'k-'); slab{end+1} = 'Spline';
  p(end+1) = plot(sval(end)/1e6,mur_ms(end,end,end),'k--'); slab{end+1} = 'Spline extrap.';
  l = legend(p, slab, 'Location', 'SouthEast'); set(l, 'FontSize', 12); legend boxoff
  xlabel('Stress {\it\sigma} (MPa)', 'FontSize', 14);
  ylabel('Relative permeability', 'FontSize', 14);
  ymin = min(min([mur_ms(:,:); mur_i(:,:); mur_i_xt(:,:)]));
  ymax = max(max([mur_ms(:,:); mur_i(:,:); mur_i_xt(:,:)]));
  ylim([ymin ymax]);
figure(2);
  xlabel('Stress {\it\sigma} (MPa)', 'FontSize', 14);
  ylabel('Magnetostriction {\it\lambda}_{xx} (ppm)', 'FontSize', 14);
  ymin = min(min([lxx_ms(:,:); lxx_i(:,:); lxx_i_xt(:,:)]))*1e6;
  ymax = max(max([lxx_ms(:,:); lxx_i(:,:); lxx_i_xt(:,:)]))*1e6;
  ylim([ymin ymax]);
figure(3);
  p2(end+1) = plot(sval(end)/1e6,lyy_ms(end,end,end)*1e6,'k-'); slab2{end+1} = 'Spline';
  p2(end+1) = plot(sval(end)/1e6,lyy_ms(end,end,end)*1e6,'k--'); slab2{end+1} = 'Spline extrap.';
  l = legend(p2, slab2); set(l, 'FontSize', 12); legend boxoff
  xlabel('Stress {\it\sigma} (MPa)', 'FontSize', 14);
  ylabel('Magnetostriction {\it\lambda}_{yy} (ppm)', 'FontSize', 14);
  ymin = min(min([lyy_ms(:,:); lyy_i(:,:); lyy_i_xt(:,:); lzz_ms(:,:); lzz_i(:,:); lzz_i_xt(:,:)]))*1e6;
  ymax = max(max([lyy_ms(:,:); lyy_i(:,:); lyy_i_xt(:,:); lzz_ms(:,:); lzz_i(:,:); lzz_i_xt(:,:)]))*1e6;
  ylim([ymin ymax]);
  annotation(gcf,'doublearrow',[0.285714285714286 0.283928571428571],...
    [0.645238095238095 0.55]);
  annotation(gcf,'doublearrow',[0.208928571428571 0.210714285714286],...
    [0.337095238095238 0.428571428571429]);
figure(4);
  xlabel('Stress {\it\sigma} (MPa)', 'FontSize', 14);
  ylabel('Magnetostriction {\it\lambda}_{zz} (ppm)', 'FontSize', 14);
  ylim([ymin ymax]);
  annotation(gcf,'doublearrow',[0.2 0.173214285714286],...
    [0.565666666666667 0.461904761904762]);
  annotation(gcf,'doublearrow',[0.292857142857143 0.319642857142857],...
    [0.435714285714286 0.495238095238095]);

