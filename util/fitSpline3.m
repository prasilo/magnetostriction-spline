function s = fitSpline3(ordr, u, v, w, psi_u, psi_v, psi_w, flag_rank, weight_u, weight_v, weight_w)
% Fit trivariate spline based on given grid points and partial derivatives

  % Calculate ranks or not (may be slow)
  if nargin < 8
    flag_rank = 0;
  end
  
  % For weighting u, v or w in the least-squares fitting
  if nargin < 9
    weight_u = 1;
  end
  if nargin < 10
    weight_v = 1;
  end
  if nargin < 11
    weight_w = 1;
  end

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
      
  % Create matrix an right-hand side for the least-squuares problem
  mat = [weight_u*umat; weight_v*vmat; weight_w*wmat];
  rhs = [weight_u*urhs; weight_v*vrhs; weight_w*wrhs];
  
  % Since only the derivatives are fitted, the value of one coefficient has
  % to be fixed -> remove one column. Get also an estimate of the ranks and
  % condition numbers before and after removing the column.
  if flag_rank
    kb = condest(mat'*mat);
    rb = rank(full(mat'));
  end
  mat(:,1) = [];
  if flag_rank
    ka = condest(mat'*mat);
    ra = rank(full(mat'));
      
    % Display ranks and condition numbers
    fprintf('Full matrix: rb = %d, kb = %.3g, one column eliminated: ra = %d, ka = %.3g. ', rb, kb, ra, ka);   
  end
  
  % Solve spline coefficients
  c = lsqlin(mat'*mat, mat'*rhs);

  % Set first coefficient to zero and create spline in pp-form
  c = [0; c];
  s.coefs = reshape(full(c), size(s.coefs));
  s = fn2fm(s, 'pp');
  
  % Separate spline for extrapolation
  s.sx = fnxtr(s,[3 3 3]);
  