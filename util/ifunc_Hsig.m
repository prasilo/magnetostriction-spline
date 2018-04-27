function [B,lam] = ifunc_Hsig(s, H, sig)
% Evaluate the spline-based invariant model with H and sigma as the state
% variables

  % Check if 2- or 3-variate spline is provided
  if length(s.order) == 2    
    flag3 = 0;
    Hscale = s.Hscale;
    sxx_scale = s.sscale;
  else
    flag3 = 1;
    Hscale = s.Hscale;
    sxx_scale = s.sxx_scale;
    sxy_scale = s.sxy_scale;
  end;

  % Input variables
  Hx = H(1,:);
  Hy = H(2,:);
  Hz = H(3,:);
  sxx = sig(1,:);
  syy = sig(2,:);
  szz = sig(3,:);
  syz = sig(4,:);
  szx = sig(5,:);
  sxy = sig(6,:);

  % Invariants and auxiliary variables
  vars = {Hx, Hy, Hz, sxx, syy, szz, syz, szx, sxy};
  u = fu(vars{:})/Hscale;
  v = fv(vars{:})/sxx_scale;
  if flag3
    w = fw(vars{:})/sxy_scale;
  end;

  % Derivatives with respect to H
  [dudhx,dudhy,dudhz] = fdudh(vars{:});
  [dvdhx,dvdhy,dvdhz] = fdvdh(vars{:});
  if flag3
    [dwdhx,dwdhy,dwdhz] = fdwdh(vars{:});
  end;
  dudhx = dudhx/Hscale;
  dudhy = dudhy/Hscale;
  dudhz = dudhz/Hscale;
  dvdhx = dvdhx/sxx_scale;
  dvdhy = dvdhy/sxx_scale;
  dvdhz = dvdhz/sxx_scale;
  if flag3
    dwdhx = dwdhx/sxy_scale;
    dwdhy = dwdhy/sxy_scale;
    dwdhz = dwdhz/sxy_scale;
  end;

  % Derivatives with respect to sigma
  [dudsxx,dudsyy,dudszz,dudsyz,dudszx,dudsxy] = fduds(vars{:});
  [dvdsxx,dvdsyy,dvdszz,dvdsyz,dvdszx,dvdsxy] = fdvds(vars{:});
  if flag3
    [dwdsxx,dwdsyy,dwdszz,dwdsyz,dwdszx,dwdsxy] = fdwds(vars{:});
  end;
  dudsxx = dudsxx/Hscale;
  dudsyy = dudsyy/Hscale;
  dudszz = dudszz/Hscale;
  dudsyz = dudsyz/Hscale;
  dudszx = dudszx/Hscale;
  dudsxy = dudsxy/Hscale;
  dvdsxx = dvdsxx/sxx_scale;
  dvdsyy = dvdsyy/sxx_scale;
  dvdszz = dvdszz/sxx_scale;
  dvdsyz = dvdsyz/sxx_scale;
  dvdszx = dvdszx/sxx_scale;
  dvdsxy = dvdsxy/sxx_scale;
  if flag3
    dwdsyz = dwdsyz/sxy_scale;
    dwdszx = dwdszx/sxy_scale;
    dwdsxy = dwdsxy/sxy_scale;

    % Correct derivatives for w = 0
    dwdhx(w == 0) = 0;
    dwdhy(w == 0) = 0;
    dwdhz(w == 0) = 0;
    dwdsxx(w == 0) = 0;
    dwdsyy(w == 0) = 0;
    dwdszz(w == 0) = 0;
    dwdsyz(w == 0) = 0;
    dwdszx(w == 0) = 0;
    dwdsxy(w == 0) = 0;
  end;

  % Evaluate spline
  if flag3 == 0
    sdu = fnder(s, [1 0]);
    sdv = fnder(s, [0 1]);
    pu = fnval(sdu, [u; v]);
    pv = fnval(sdv, [u; v]);

    % Calculate B
    Bx = (pu.*dudhx + pv.*dvdhx);
    By = (pu.*dudhy + pv.*dvdhy);
    Bz = (pu.*dudhz + pv.*dvdhz);

    % Calculate lambda
    lxx = (pu.*dudsxx + pv.*dvdsxx);
    lyy = (pu.*dudsyy + pv.*dvdsyy);
    lzz = (pu.*dudszz + pv.*dvdszz);
    lyz = (pu.*dudsyz + pv.*dvdsyz);
    lzx = (pu.*dudszx + pv.*dvdszx);
    lxy = (pu.*dudsxy + pv.*dvdsxy);

  else
    sdu = fnder(s, [1 0 0]);
    sdv = fnder(s, [0 1 0]);
    sdw = fnder(s, [0 0 1]);
    pu = fnval(sdu, [u; v; w]);
    pv = fnval(sdv, [u; v; w]);
    pw = fnval(sdw, [u; v; w]);  

    % Calculate B
    Bx = (pu.*dudhx + pv.*dvdhx + pw.*dwdhx);
    By = (pu.*dudhy + pv.*dvdhy + pw.*dwdhy);
    Bz = (pu.*dudhz + pv.*dvdhz + pw.*dwdhz);

    % Calculate lambda
    lxx = (pu.*dudsxx + pv.*dvdsxx + pw.*dwdsxx);
    lyy = (pu.*dudsyy + pv.*dvdsyy + pw.*dwdsyy);
    lzz = (pu.*dudszz + pv.*dvdszz + pw.*dwdszz);
    lyz = (pu.*dudsyz + pv.*dvdsyz + pw.*dwdsyz);
    lzx = (pu.*dudszx + pv.*dvdszx + pw.*dwdszx);
    lxy = (pu.*dudsxy + pv.*dvdsxy + pw.*dwdsxy);  
  end;

  % Output vectors
  B = [Bx; By; Bz];
  lam = [lxx; lyy; lzz; lyz; lzx; lxy];

  % Extrapolated values to nan
  iu = (u < min(s.breaks{1})) | (u > max(s.breaks{1})); 
  iv = (v < min(s.breaks{2})) | (v > max(s.breaks{2})); 
  if flag3
    iw = (w < min(s.breaks{3})) | (w > max(s.breaks{3})); 
    ixtr = (iu | iv | iw);
  else
    ixtr = (iu | iv);
  end;
  B(:,ixtr) = nan;
  lam(:,ixtr) = nan;
