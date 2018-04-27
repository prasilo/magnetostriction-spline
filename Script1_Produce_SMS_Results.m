% Produces reference data with the simplified multiscale model.

clear all
close all
clc
addpath util

% Bi-variate data
Hx = linspace(0,6000,40);
sigxx = linspace(-100,100,11)*1e6;
col = jet(length(sigxx));
ns = length(sigxx);
nh = length(Hx);
for is = 1:ns
  sig = [sigxx(is) 0 0; 0 0 0; 0 0 0];
  for ih = 1:nh
    fprintf('%d\\%d %d\\%d\n', is, ns, ih, nh);

    % Run simplified multiscale model
    H = [Hx(ih); 0; 0];
    [B,lam] = smsfunc(H, sig);

    % Values to tables
    Bx(ih,is) = B(1);
    lamxx(ih,is) = lam(1,1);
    lamyy(ih,is) = lam(2,2);
    lamzz(ih,is) = lam(3,3);
  end;

  % Plot
  figure(10);
  subplot(2,1,1);
    hold on;
    plot(Hx, Bx(:,is), '.-', 'Color', col(is,:));
  subplot(2,1,2);
    plot(Bx(:,is), lamxx(:,is)-lamxx(1,is),'.-', 'Color', col(is,:))
    hold on
  drawnow;
end;
save .\data\sms_data2 Hx sigxx Bx lamxx lamyy lamzz
  
% Trivariate data
clear all;
Hx = linspace(0,6000,40);
sigxx = linspace(-100,100,11)*1e6;
sigxy = linspace(-100,100,11)*1e6;
col = jet(length(sigxx));
nsxx = length(sigxx);
nsxy = length(sigxy);
nh = length(Hx);
for isxx = 1:nsxx
  for isxy = 1:nsxy
    sig = [sigxx(isxx) sigxy(isxy) 0; sigxy(isxy) 0 0; 0 0 0];
    for ih = 1:nh
      fprintf('%d\\%d %d\\%d %d\\%d\n', isxx, nsxx, isxy, nsxy, ih, nh);
      H = [Hx(ih); 0; 0];

      % Run simplified multiscale model
      [B,lam] = smsfunc(H, sig);

      % Values to tables
      Bx(ih,isxx,isxy) = B(1);
      lamxx(ih,isxx,isxy) = lam(1,1);
      lamxy(ih,isxx,isxy) = lam(1,2);
      lamyy(ih,isxx,isxy) = lam(2,2);
      lamzz(ih,isxx,isxy) = lam(3,3);
    end;

    % Plot
    figure(100+isxx);
    subplot(3,1,1);
      hold on;
      plot(Hx, Bx(:,isxx,isxy), '.-', 'Color', col(isxy,:));
    subplot(3,1,2);
      plot(Bx(:,isxx,isxy), lamxx(:,isxx,isxy)-lamxx(1,isxx,isxy),'.-', 'Color', col(isxy,:))
      hold on
    subplot(3,1,3);
      plot(Bx(:,isxx,isxy), lamxy(:,isxx,isxy)-lamxy(1,isxx,isxy),'.-', 'Color', col(isxy,:))
      hold on
    drawnow;
  end;
end;
save .\data\sms_data3 Hx sigxx sigxy Bx lamxx lamxy lamyy lamzz