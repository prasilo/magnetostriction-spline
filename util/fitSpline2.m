function s = fitSpline2(ordr, u, v, psi_u, psi_v, weight_u, weight_v, uknots, vknots)
% Fit bivariate spline based on given grid points and partial derivatives

  % For weighting u or v in the least-squares fitting
  if nargin < 6
    weight_u = 1;
  end;
  if nargin < 7
    weight_v = 1;
  end;
  if nargin < 8
    uknots = optknt(u,ordr);
  end;
  if nargin < 9
    vknots = optknt(v,ordr);
  end;

  % Number of points
  Nu = length(u);
  Nv = length(v);

  % Spline collocation matrices, Equation (10)
  Au = spcol(uknots, ordr, u);
  Av = spcol(vknots, ordr, v);

  % Collocation matrices for the partial derivatives:
  % These map the partial derivative spline coefficients into the
  % partial derivative values at the measurement points
  Adu = spcol(uknots(2:end-1), ordr-1, u);
  Adv = spcol(vknots(2:end-1), ordr-1, v);

  % Coefficient differentiation matrices:
  % These map the energy-density spline coefficients into the partial
  % derivative spline coefficients
  Dcu = full(DerivBKnotDeriv(uknots, ordr, 1, ones(length(uknots), 1)));
  Dcv = full(DerivBKnotDeriv(vknots, ordr, 1, ones(length(vknots), 1)));

  % Differentation matrices, Equation (11):
  % These map the energy-density spline coefficients into the partial
  % derivative values at the measurement points
  Du = Adu*Dcu;
  Dv = Adv*Dcv;
  
  % Make a dummy spline structure for these knots
  % The coefficients will later be replaced the identified ones
  s = spap2({uknots,vknots}, ordr, {u,v}, 0*psi_u); 
  
%%% Assemble system and solve coefficients

  % Size of coefficient matrix
  cu = length(uknots)-ordr;

  % Summing indices calculated beforehand to avoid for-loops
  l  = 1 : size(Au,2);
  m  = 1 : size(Av,2);
  [l,m] = ndgrid(l,m); l = l(:); m = m(:);

  % i, j and k indices and row index corresponding to the equation index
  [i,j] = ndgrid(1:Nu,1:Nv); i = i(:); j = j(:);
  row = max(i)*(j-1) + i;

  % u matrix and rhs
  indu = (m-1)*cu + l;
  umat(row,indu) = Av(j,m).*Du(i,l);
  umat = sparse(umat);
  urhs(row,1) = sparse(psi_u(row));

  % v matrix and rhs
  indv = (m-1)*cu + l;
  vmat(row,indv) = Au(i,l).*Dv(j,m);
  vmat = sparse(vmat);
  vrhs(row,1) = sparse(psi_v(row));

  % Instead of solving A(i,:)*x = bi, solve A(i,:)/bi*x = 1,
  % so that all points have the same weight (except if bi = 0)
  oo = ones(1,size(umat,2));
  sc = 0*urhs+1; ind = (abs(urhs) > 1e-6); sc(ind,1) = urhs(ind); urhs(ind) = 1; sc = sc*oo; umat = umat./sc;
  sc = 0*vrhs+1; ind = (abs(vrhs) > 1e-6); sc(ind,1) = vrhs(ind); vrhs(ind) = 1; sc = sc*oo; vmat = vmat./sc;
  mat = [weight_u*umat; weight_v*vmat];
  rhs = [weight_u*urhs; weight_v*vrhs];
  
  % Since only the derivatives are fitted, the value of one coefficient has
  % to be fixed
  mat(1,:) = 0; mat(1,1) = 1; rhs(1) = 1;

  % Solve spline coefficients and change to pp-form
  s.coefs = reshape(full(lsqlin(mat,rhs)), size(s.coefs));
  s = fn2fm(s, 'pp');
  
  % Extrapolate (separate spline) (3 means that the energy is extrapolated
  % as a second-order polynomial), change to pp-form
  s.sx = fnxtr(s,[3 3]);
  s.sx = fn2fm(s.sx, 'pp');
