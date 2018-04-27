function s = fitSpline3(ordr, u, v, w, psi_u, psi_v, psi_w, weight_u, weight_v, weight_w)
% Fit trivariate spline based on given grid points and partial derivatives

  % For weighting u, v or w in the least-squares fitting
  if nargin < 8
    weight_u = 1;
  end;
  if nargin < 9
    weight_v = 1;
  end;
  if nargin < 10
    weight_w = 1;
  end;

  % Number of points
  Nu = length(u);
  Nv = length(v);
  Nw = length(w);
  
  % Knots
  uknots = optknt(u, ordr);
  vknots = optknt(v, ordr);
  wknots = optknt(w, ordr);

  % Spline collocation matrices, Equation (10)
  Au = spcol(uknots, ordr, u);
  Av = spcol(vknots, ordr, v);
  Aw = spcol(wknots, ordr, w);  
  
  % Collocation matrices for the partial derivatives:
  % These map the partial derivative spline coefficients into the
  % partial derivative values at the measurement points
  Adu = spcol(uknots(2:end-1), ordr-1, u);
  Adv = spcol(vknots(2:end-1), ordr-1, v);
  Adw = spcol(wknots(2:end-1), ordr-1, w);

  % Coefficient differentiation matrices:
  % These map the energy-density spline coefficients into the partial
  % derivative spline coefficients
  Dcu = DerivBKnotDeriv(uknots, ordr, 1, ones(length(uknots), 1));
  Dcv = DerivBKnotDeriv(vknots, ordr, 1, ones(length(vknots), 1));
  Dcw = DerivBKnotDeriv(wknots, ordr, 1, ones(length(wknots), 1));

  % Differentation matrices, Equation (11):
  % These map the energy-density spline coefficients into the partial
  % derivative values at the measurement points
  Du = Adu*Dcu;
  Dv = Adv*Dcv;
  Dw = Adw*Dcw;

  % Make a dummy spline structure for these knots
  % The coefficients will later be replaced by the identified ones
  s = spap2({uknots,vknots,wknots}, ordr, {u,v,w}, 0*psi_u);  
  
%%% Assemble system and solve coefficients
  
  % Size of coefficient matrix
  cu = length(uknots)-ordr;
  cv = length(vknots)-ordr;
  cuv = cu*cv;

  % Summing indices calculated beforehand to avoid for-loops
  l  = 1 : size(Au,2);
  m  = 1 : size(Av,2);
  n  = 1 : size(Aw,2);
  [l,m,n] = ndgrid(l,m,n); l = l(:); m = m(:); n = n(:);

  % i, j and k indices and row index corresponding to the equation index
  [i,j,k] = ndgrid(1:Nu,1:Nv,1:Nw); i = i(:); j = j(:); k = k(:);
  row = max(j)*max(i)*(k-1) + max(i)*(j-1) + i;

  % u matrix and rhs
  indu = (n-1)*cuv + (m-1)*cu + l;
  umat(row,indu) = Av(j,m).*Aw(k,n).*Du(i,l);
  umat = sparse(umat);
  urhs(row,1) = sparse(psi_u(row));

  % v matrix and rhs
  indv = (n-1)*cuv + (m-1)*cu + l;
  vmat(row,indv) = Au(i,l).*Aw(k,n).*Dv(j,m);
  vmat = sparse(vmat);
  vrhs(row,1) = sparse(psi_v(row));

  % z matrix and rhs
  indw = (n-1)*cuv + (m-1)*cu + l;
  wmat(row,indw) = Au(i,l).*Av(j,m).*Dw(k,n);
  wmat = sparse(wmat);
  wrhs(row,1) = sparse(psi_w(row));
      
  % Instead of solving A(i,:)*x = bi, solve A(i,:)/bi*x = 1,
  % so that all points have the same weight (except if bi = 0)
  oo = ones(1,size(umat,2));
  sc = 0*urhs+1; ind = (abs(urhs) > 1e-6); sc(ind,1) = urhs(ind); urhs(ind) = 1; sc = sc*oo; umat = umat./sc;
  sc = 0*vrhs+1; ind = (abs(vrhs) > 1e-6); sc(ind,1) = vrhs(ind); vrhs(ind) = 1; sc = sc*oo; vmat = vmat./sc;
  sc = 0*wrhs+1; ind = (abs(wrhs) > 1e-6); sc(ind,1) = wrhs(ind); wrhs(ind) = 1; sc = sc*oo; wmat = wmat./sc;
  mat = [weight_u*umat; weight_v*vmat; weight_w*wmat];
  rhs = [weight_u*urhs; weight_v*vrhs; weight_w*wrhs];
  
  % Since only the derivatives are fitted, the value of one coefficient has
  % to be fixed
  mat(1,:) = 0; mat(1,1) = 1; rhs(1) = 1;

  % Solve spline coefficients and change to pp-form
  s.coefs = reshape(full(lsqlin(mat,rhs)), size(s.coefs));
  s = fn2fm(s, 'pp');

  % Extrapolate (separate spline) (3 means that the energy is extrapolated
  % as a second-order polynomial), change to pp-form
  s.sx = fnxtr(s,[3 3 3]);
  s.sx = fn2fm(s.sx, 'pp');
