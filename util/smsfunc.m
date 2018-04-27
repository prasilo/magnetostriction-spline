function [B,lam] = bfunc(H, sig)
% [1]	L. Daniel, O. Hubert, M. Rekik, A Simplified 3D Constitutive Law for Magneto-Mechanical Behavior, IEEE Trans. Magn. 51(3) (2015), 7300704.

  % From Voigt to tensor
  if min(size(sig)) == 1
    sig = sig([1 6 5; 6 2 4; 5 4 3]);
  end;

  % Parameters
  mu0 = 4*pi*1e-7;
  Ms = 1.45e6;
  lams = 2/3*10e-6;
  As = 1.8e-3;
  chi0 = mu0*As*Ms^2/3;
  
  % Local magnetization, Equation (5)
  f_Malfa1= @(theta,phi) Ms*sin(phi).*cos(theta); 
  f_Malfa2= @(theta,phi) Ms*sin(phi).*sin(theta);
  f_Malfa3= @(theta,phi) Ms*cos(phi);

  % Local magnetostriction, Equation (6)
  f_lalfa1  = @(theta,phi)   3/2*lams*((sin(phi).*cos(theta)).^2-1/3);
  f_lalfa2  = @(theta,phi)   3/2*lams*((sin(phi).*cos(theta)).*(sin(phi).*sin(theta)));
  f_lalfa3  = @(theta,phi)   3/2*lams*((sin(phi).*cos(theta)).*(cos(phi)));
  f_lalfa4  = @(theta,phi)   3/2*lams*((sin(phi).*cos(theta)).*(sin(phi).*sin(theta)));
  f_lalfa5  = @(theta,phi)   3/2*lams*((sin(phi).*sin(theta)).^2-1/3) ;
  f_lalfa6  = @(theta,phi)   3/2*lams*((sin(phi).*sin(theta)).*(cos(phi)));
  f_lalfa7  = @(theta,phi)   3/2*lams*((sin(phi).*cos(theta)).*(cos(phi))) ;
  f_lalfa8  = @(theta,phi)   3/2*lams*((sin(phi).*sin(theta)).*cos(phi))  ;
  f_lalfa9  = @(theta,phi)   3/2*lams*((cos(phi)).^2-1/3);

  % Magnetic part of the energy, Equation (2)
  f_Wmag   = @(theta,phi) -mu0*(H(1)*f_Malfa1(theta,phi)+H(2)*f_Malfa2(theta,phi)+H(3)*f_Malfa3(theta,phi));  % equ (2)

  % Mechanical part of the energy, Equation (3)
  f_Wmek   = @(theta,phi) -(sig(1,1)*f_lalfa1(theta,phi)+sig(1,2)*f_lalfa2(theta,phi)+sig(1,3)*f_lalfa3(theta,phi) ...
                           +sig(2,1)*f_lalfa4(theta,phi)+sig(2,2)*f_lalfa5(theta,phi)+sig(2,3)*f_lalfa6(theta,phi) ...
                           +sig(3,1)*f_lalfa7(theta,phi)+sig(3,2)*f_lalfa8(theta,phi)+sig(3,3)*f_lalfa9(theta,phi)); % equ (3)

  % Total energy, Equation (1)
  f_Walfa  = @(theta,phi)  f_Wmag(theta,phi) + f_Wmek(theta,phi);

  % Probability, this already requires numerical integration, Equation (7)
  funcint = @(theta,phi) sin(phi).*exp(-As*f_Walfa(theta,phi));
  intfuncin = quad2d(funcint,0,2*pi,0,pi, 'abstol', 1e-12, 'reltol', 1e-6);
  falfa = @(theta,phi) exp(-As*f_Walfa(theta,phi))./ intfuncin;

  % Integrands for magnetization
  f1=@(theta,phi) sin(phi).*f_Malfa1(theta,phi).*falfa(theta,phi);
  f2=@(theta,phi) sin(phi).*f_Malfa2(theta,phi).*falfa(theta,phi);
  f3=@(theta,phi) sin(phi).*f_Malfa3(theta,phi).*falfa(theta,phi);

  % Integrate for magnetization, Equation (9)
  Mx =  quad2d(f1, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  My =  quad2d(f2, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6); 
  Mz =  quad2d(f3, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6); 
  M = [Mx; My; Mz];

  % Integrands for magnetostriction
  fe1=@(theta,phi) sin(phi).*f_lalfa1(theta,phi).*falfa(theta,phi);
  fe2=@(theta,phi) sin(phi).*f_lalfa2(theta,phi).*falfa(theta,phi);
  fe3=@(theta,phi) sin(phi).*f_lalfa3(theta,phi).*falfa(theta,phi);    
  fe5=@(theta,phi) sin(phi).*f_lalfa5(theta,phi).*falfa(theta,phi);
  fe6=@(theta,phi) sin(phi).*f_lalfa6(theta,phi).*falfa(theta,phi);    
  fe9=@(theta,phi) sin(phi).*f_lalfa9(theta,phi).*falfa(theta,phi);            

  % Integrate for magnetostrictions, Equation (10)
  lam1 =  quad2d(fe1, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam2 =  quad2d(fe2, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam3 =  quad2d(fe3, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam4 =  lam2;   
  lam5 =  quad2d(fe5, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam6 =  quad2d(fe6, 0, 2*pi, 0, pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam7 =  lam3;   
  lam8 =  lam6;   
  lam9 =  quad2d(fe9, 0, 2*pi,0,pi, 'abstol', 1e-12, 'reltol', 1e-6);   
  lam=[lam1 lam2 lam3; lam4 lam5 lam6; lam7 lam8 lam9];       

  % Flux density
  B = mu0*(H + M);
