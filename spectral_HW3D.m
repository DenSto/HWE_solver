function spectral_HW3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script to simulate the Hasegawa-Wakatani equation in two and           %
% three dimensions. Includes phase-shift-periodic boundary conditions           %
% to model an irrational flux surface, needed to kill convective cells          %
% in the three dimensional case.                                                %
%                                                                               %
% Based off an MIT code originally made by Jean-Christophe Nave.                %
% Modified by Denis St-Onge                                                     %
%                                                                               %
% Laplacian(phi) = w                                                            %
% u = phi_y                                                                     %
% v =-phi_x                                                                     %
%                                                                               %
% Colorbrewer and Ander Biguri's Perceptually Uniform Colormaps are used        %
%                                                                               %
% Script to write to gnuplot binaries can be obtained here                      %
% http://www.gnuplotting.org/code/save_binary_matrix.m                          %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dto lg M r Q1 Q2 f1 f2 f3 isreal; 
clearvars ;
clear ETDRK4; clear ETDRK3; clear ETDRK2;
basename='HWE_2d/';


cm_redblue = flip(cbrewer('div', 'RdBu',129));
cm_inferno = inferno();
c_map      = cm_inferno;
c_maprb    = cm_redblue;

mkdir(basename);
mkdir([basename 'plots']); mkdir([basename 'data']);
delete([basename 'log.txt']);

scale = 1;                 % quick scale of the linear terms
mu   = scale*[0.0];        % friction, hyperfriction (hypoviscosity), etc...
nu   = scale*[0,0.01];  % viscosity, hyperviscosity, 6th order hyperdisssipation, etc...
chi = 0.125;            % adiabaticity parameter chi (alpha in 2D or NZ_real = 1)
kappa = 1;              % density gradient paramger

LX=2*pi*10;             % X length
LY=2*pi*10;             % Y length
LZ=2*pi;                % Z scale
NX_real = 128;          % resolution in x before padding
NY_real = 128;          % resolution in y before padding
NZ_real = 32;           % resolution in z before padding (NZ_real = 1 for 2D)
irrational = true;      % If true, apply phase shift across parallel (Z) direction
coprime = 37;           % Phase shift factor is coprime / NY. If negative, a random coprime of 
                        % NY is chosen.
 
dt = 1e-2;              % time step. Should start small as CFL updated can pick up the pace
pert_size = 1e-3;       % size of perturbation

TF=15000.0;             % final time
iF=100000000;           % final iteration, whichever occurs first
iRST=200000;            % cadence of writing restart dumps

initial_condition='random';   %'random' or 'restart' 

i_report=100;           % cadence of writing to log.txt
en_print=10;            % cadence of writing to the energy file
AB_order=3;             % Adams-Bashforth order (1-4)
TSCREEN=1000;           % sreen update interval time (NOTE: plotting is usually slow)
simulation_type='NL';   % NL or L (nonlinear, quasilinear, flux-balanced and linear respectively)
padding = true;         % 3/2 padding,spspec otherwise 2/3 truncation (latter doesn't work yet)
with_plotting = true;   % make plots
save_plots = true;      % save plots to file

cfl_cadence=1;          % cadence of recalculating timestep
cfl=0.5;                % CFL number
safety=0.8;             % CFL safety number
max_dt=1e-1;            % maximum timestep

%rng(707296708);
rng('shuffle');
s=rng;

% print log file.
diary([basename 'log.txt']) 
diary on;

fig1=0;
if(with_plotting)
  fig1=figure(1);
end


NX = NX_real;
NY = NY_real;
NZ = NZ_real;
if(padding)
  NX = 3 * (NX + mod(NX, 2)) / 2;
  NY = 3 * (NY + mod(NY, 2)) / 2;
  if(NZ > 1) 
    NZ = 3 * (NZ - mod(NZ, 2)) / 2;
  end
end

if(coprime < 0)
   coprime = get_coprime;
end

%ensure parameters get printed  to log.
fprintf('mu: ');
for i = 1:length(mu), fprintf('  %.05e',mu(i)); end; fprintf('\n');
fprintf('nu: ');
for i = 1:length(nu), fprintf('  %.05e',nu(i)); end; fprintf('\n');
fprintf('chi:%.05e kappa: %.05e  irrational:%d coprime:%d\n',chi,kappa,irrational,coprime);
fprintf('LX:%.02f LY:%.02f LZ:%.02f NX:%d NY:%d Nz:%d\n',LX, LY, LZ, NX_real, NY_real, NZ_real);
fprintf('scale:%d Tf:%.01f iF:%d\n', scale, TF, iF);
fprintf('Nonlinear:%s padding:%d\n',simulation_type, padding);
fprintf('random seed:%d AB order:%d CFL step:%d\n',s.Seed, AB_order, cfl_cadence);
fprintf('safety:%f perturbation size:%.05e\n',safety, pert_size);

if(strcmp(initial_condition,'restart'))
  energyFile = fopen([basename 'energy.dat'],'a');
else
  energyFile = fopen([basename 'energy.dat'],'w');
end

fprintf(energyFile,'# [1] t  [2] dt  [3] energy  [4] enstrophy  [5] ZF    [6] DW    [7] flux\n');


dx = LX/NX;
dy = LY/NY;
dz = LZ/NZ;
dkx = 2*pi/LX;
dky = 2*pi/LY;
dkz = 2*pi/LZ;
LXnum = 0:dx:LX;
LYnum = 0:dy:LY;
if(NZ == 1)
  z = 0;
else
  z = repmat(reshape((1:NZ)*dz,1,1,NZ),NX,NY,1);
end

minkx = -(NX/2 - 0)*dkx;
maxkx =  (NX/2 - 1)*dkx;
minky = -(NY/2 - 0)*dky;
maxky =  (NY/2 - 1)*dky;
minkz = -(NZ/2 - 0)*dkz;
maxkz =  (NZ/2 - 1)*dkz;
kxnum = minkx:dkx:maxkx;
kynum = minky:dky:maxky;

t=0.;
i=0;

cp1=0;
cp2=0;
cp3=0;
cn1=0;
cn2=0;
cn3=0;
dt1=0;
dt2=0;
dt3=0;

enk=0;
rsk=0;
rek=0;
trk=0;
phi_hat=0; %this is our field.
I=sqrt(-1);

[kx,ky,kz] = ndgrid(ifftshift(minkx:dkx:maxkx), ...
                    ifftshift(minky:dky:maxky), ...
                    ifftshift(minkz:dkz:maxkz));

if(NZ == 1)
  if(irrational)
    kz = ones(NX,NY,NZ);
    kz(:,1,1) = 0;
  else
    kz = 1;
  end
elseif(irrational)
  kz = kz + ky.*(coprime*dy/LZ);  
end
                  
if(irrational)
   zfac = exp(-I.*ky.*z.*(coprime*dy/LZ));
  izfac = 1.0./zfac;
else
  zfac  = 1.0;
  izfac = 1.0;
end

% Cutting of frequencies using the 2/3 rule. 
% Discard asymmetric N/2 term. Otherwise reality is not enforced.
if(NZ==1)
  dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY;
else
  dealias=abs(kx/dkx) <1/3*NX & abs(ky/dky) <1/3*NY & abs(kz/dkz) <1/3*NZ;
end
dealias(1,1,1)=0;

ksquare = kx.^2 + ky.^2;                % Laplacian in Fourier space 
iksquare = 1.0./ksquare;
iksquare(1,1,:) = 0.0;
kmu = build_friction(mu);               % Friction 
knu = build_viscosity(nu);              % Viscosity


T11 = -ksquare;
T12 = 0;
T21 = 0;
T22 = 1.0;
L11 = -chi.*(kz.^2) - (knu + kmu).*(-ksquare);
L12 =  chi.*(kz.^2);
L21 = -kappa*I*ky - chi.*(kz.^2);
L22 = - (knu + kmu)  + chi.*(kz.^2);

B11=0; B12=0; B21=0; B22=0;
C11=0; C12=0; C21=0; C22=0;


n_hat   = zeros(NX,NY,NZ);
phi_hat = zeros(NX,NY,NZ);

% Define initial vorticity distribution
switch lower(initial_condition)
    case {'random'}
      phi=pert_size*randn(NX,NY,NZ);%normally 5e-2
      ne =pert_size*randn(NX,NY,NZ);%normally 5e-2

      phi_hat=fftz(fft2(phi))./(ksquare+kz.^2+0.0000001);
      n_hat =fftz(fft2(ne));
    case {'restart'}
      fileID=fopen('final_state.bin','r');
      t=fread(fileID,1,'double');
      dt=fread(fileID,1,'double');
      i=fread(fileID,1,'int');
      enk=fread(fileID,1,'int');
      trk=fread(fileID,1,'int');
      rsk=fread(fileID,1,'int');
      rek=fread(fileID,1,'int');
      phi_hat=fread(fileID,[NX NY NZ],'double')+I*fread(fileID,[NX NY NZ],'double');  
      n_hat  =fread(fileID,[NX NY NZ],'double')+I*fread(fileID,[NX NY NZ],'double');            
     
      fclose(fileID);        
    otherwise
      disp('Unknown initial conditions');
      return
end

if(padding) % ensure there's no energy in the padded region of k-space.
  phi_hat=dealias.*phi_hat;
  n_hat  =dealias.*n_hat;
end

phi_hat= enforceReality(phi_hat);
n_hat  = enforceReality(n_hat);

u=0;
v=0;

calc_matrix;

tic
while t<TF && i<iF
    if(any(isnan(phi_hat(:)))) % Something really bad has happened.
        fprintf('Divergence at iteration: %d\n',i);
        return;
    end
    if(mod(i,i_report)==0)
       fprintf('iteration: %d    dt:%.02e     t:%.03e     step time:%.1d s\n',i,dt,t,toc); 
       tic;
       diary off; %flush diary
       diary on;
    end
    if (mod(i,en_print)== 0) 
        outputEnergy();
    end
    
    
   if (mod(i,TSCREEN)== 0 && with_plotting) 
        plotfunc();
    end
    
    if (mod(i,iRST) == 0)
      dump(sprintf('%srestart.bin',basename));
      rek=rek+1;
    end
    
    [conv_p_hat,conv_n_hat] = nonlinear(phi_hat,n_hat,simulation_type); 
 
    if(mod(i+1,cfl_cadence)==0 && i > 3) % compute new timestep from CFL condition.
      abs_u=abs(u);
      abs_v=abs(v);
      maxV= max(abs_u(:))/dx + max(abs_v(:))/dy;
      if(maxV>0)
          new_dt=1/maxV;
      else
          new_dt=inf;
      end
      target_dt=min(cfl*new_dt,min(4.0*dt/safety,max_dt/safety));
      if(target_dt < dt)
        disp('WARNING: New dt fell below safety.')
      end
      if target_dt < dt/safety || target_dt > 3.2*dt
        dt = max_dt/safety;
        while dt > target_dt 
          dt=dt/2.0;
        end
        dt=safety*dt;
        calc_matrix;
      end
    end    
    
    AB1=1.0; AB2=0; AB3=0; AB4=0;

    if (i < 1 || AB_order == 1) %Forward-Euler to generate history. Should run a small time-step at first. 
    elseif (i < 2 || AB_order == 2)
      w1=dt1/dt;
      AB1=(1.0 + 0.5*w1);
      AB2=0.5*w1;
    elseif (i < 3 || AB_order == 3)
      w1=dt1/dt;
      w2=dt2/dt;
      AB1=(2.0 + 3.0*w2 + 6.0*w1*(1+w1+w2))/(6.0*w1*(w1+w2));
      AB2=(2.0 + 3.0*w1 + 3.0*w2)/(6.0*w1*w2);
      AB3=(2.0 + 3.0*w1)/(6.0*w2*(w1+w2));
    elseif (i < 4 || AB_order == 4)
      w1=dt1/dt;
      w2=dt2/dt;
      w3=dt3/dt;
      AB1=(3.0 + 8.0*w2 + 4.0*w3 + 6.0*(2*w1^3 + w2*(w2+w3) +2.0*w1*(1.0+w2)*(1.0+w2+w3) + w1^2*(3.0 + 4.0*w2+2*w3)));
      AB1=AB1/(12.0*w1*(w1+w2)*(w1+w2+w3));
      AB2=(3.0 + 6.0*(w1+w2)*(w1+w2+w3) + 4.0*(2.0*(w1+w2)+w3))/(12.0*w1*w2*(w2+w3));
      AB3=(3.0 + 6.0*w1*(w1+w2+w3)+4.0*(2.0*w1 + w2+w3))/(12.0*w2*(w1+w2)*w3);
      AB4=(3.0 + 6.0*w1*(w1+w2)+4.0*(2.0*w1 + w2))/(12.0*w3*(w2+w3)*(w1+w2+w3));
    end

    Qp = C11.*phi_hat + C12.*n_hat - dt*(AB1*conv_p_hat - AB2*cp1 + AB3*cp2 - AB4*cp3);
    Qn = C21.*phi_hat + C22.*n_hat - dt*(AB1*conv_n_hat - AB2*cn1 + AB3*cn2 - AB4*cn3);

    % Compute Solution at the next step
    % Implicitly solve the linear term with 2nd order Crank-Nicholson
    phi_hat_new = B11.*Qp + B12.*Qn;
    n_hat_new   = B21.*Qp + B22.*Qn;
      
    t=t+dt;
    i=i+1;

    cp3=cp2;
    cp2=cp1;
    cp1=conv_p_hat;

    cn3=cn2;
    cn2=cn1;
    cn1=conv_n_hat;

    dt3=dt2;
    dt2=dt1;
    dt1=dt;
    
    phi_hat = dealias.*phi_hat_new;
    n_hat   = dealias.*n_hat_new;
    
    if(mod(i,100)==0)
      phi_hat = enforceReality(phi_hat);
      n_hat   = enforceReality(n_hat);
    end
end
fclose(energyFile);

%write final state 
fprintf('Simulation ended at time %.03e and iteration %d\n',t,i);
dump('final_state.bin');

% SIMULATION ENDS HERE

%%----------------------------------%%
%% Helper Function Definitions      %% 
%%----------------------------------%%

    function y=zonal_part(x)
      y = squeeze(x(:,1,1));
    end

    function v=build_friction(amp)
        v = zeros(NX,NY);
        for i1=1:length(amp)
           v = v - amp(i1)*ones(NX,NY) .*iksquare.^(i1-1); 
        end
    end
  
    function v=build_viscosity(amp)
        v = zeros(NX,NY);
        for i1=1:length(amp)
           v = v - amp(i1)*ones(NX,NY) .* (kx.^2+ky.^2).^(i1); 
        end
    end
  
    function y=enforceReality(x) % x and y are in Fourier space
        y = fftz(fft2(real(ifft2(ifftz(x)))));
    end
  
  
    function calc_matrix
  % create matrix components
      A11 = T11 + 0.5*dt*L11;
      A12 = T12 + 0.5*dt*L12;
      A21 = T21 + 0.5*dt*L21;
      A22 = T22 + 0.5*dt*L22;
      C11 = T11 - 0.5*dt*L11;
      C12 = T12 - 0.5*dt*L12;
      C21 = T21 - 0.5*dt*L21;
      C22 = T22 - 0.5*dt*L22;
      D= A11.*A22 - A12.*A21;
      B11= A22./D;
      B12=-A12./D;
      B21=-A21./D;
      B22= A11./D;
      
      B11(1,1,:) = 0.0; B12(1,1,:) = 0.0;
      B21(1,1,:) = 0.0; B22(1,1,:) = 0.0;
    end
  
   function y=get_coprime
     all_n = 1:(NY-1);
     cop_arr = (gcd(all_n,NY) == 1);
     all_co = all_n(cop_arr);
     nco = length(all_co);
     y = all_co(randi(nco));
   end
    
   function y = fftz(x)
     y =  fft(zfac.*x,[],3);
   end

   function y = ifftz(x)
     y = izfac.*ifft(x,[],3);
   end
      
    function [phi_nl,n_nl] = nonlinear(phi_h,n_h,type)
      cp_hat=0;
      cn_hat=0;
      switch upper(type)
         case {'NL'} 
              uhat = -I*ky.*phi_h;
              vhat =  I*kx.*phi_h;
              
              vor_xhat =  I*kx.*phi_h.*(-ksquare);
              vor_yhat =  I*ky.*phi_h.*(-ksquare);
              n_xhat =   I*kx.*n_h;
              n_yhat =   I*ky.*n_h;
              
              % dealiasing here truncates if not padded, otherwise it has no effect
              u  =  real(ifft2(ifftz(dealias.*uhat)));         % Compute -y derivative of stream function ==> u
              v  =  real(ifft2(ifftz(dealias.*vhat)));         % Compute  x derivative of stream function ==> v

              n_x = real(ifft2(ifftz(dealias.*n_xhat)));       % Compute  x derivative of density
              n_y = real(ifft2(ifftz(dealias.*n_yhat)));       % Compute  y derivative of density
              vor_x=real(ifft2(ifftz(dealias.*vor_xhat)));     % Compute  x derivative of vorticity
              vor_y=real(ifft2(ifftz(dealias.*vor_yhat)));     % Compute  y derivative of vorticity
              
              conv_p     = u.*vor_x + v.*vor_y;       % evaluate the convective derivative (u,v).grad(w)   
              conv_n     = u.*n_x   + v.*n_y;         % evaluate the convective derivative (u,v).grad(w)   
              
              cp_hat = fftz(fft2(conv_p));
              cn_hat = fftz(fft2(conv_n));
          case {'L'}%full linear 
              %do nothing
            otherwise
              disp('Unknown simulation type.');
              return 
      end
      if(padding)
        cp_hat    = dealias.*cp_hat;
        cn_hat    = dealias.*cn_hat;
      end
      phi_nl = cp_hat;
      n_nl   = cn_hat;
    end

    function outputEnergy()
      diary off; %flush diary
      diary on;


      w_curr  = -ksquare.*phi_hat;
      phi_curr= phi_hat;
      n_curr = n_hat;
      phi_y   = I*phi_curr.*ky;

      enstrophy = 0.5*(w_curr-n_curr).*conj(w_curr-n_curr)/(NX*NY*NZ)^2;
      energy    = 0.5*real(-conj(phi_curr).*w_curr ... 
                           +conj(n_curr).*n_curr)/(NX*NY*NZ)^2;
      flux = 0.5*real(conj(phi_y).*n_curr)/(NX*NY*NZ)^2;
     
      enstrophy_tot = sum(enstrophy(:));
      energy_tot    = sum(energy(:)); 
      ZF_energy = zonal_part(energy);
      ZF_energy_tot = sum(ZF_energy(:));
      DW_energy = energy;
      DW_energy(:,1,1) = DW_energy(:,1,1) - ZF_energy;
      DW_energy_tot = sum(DW_energy(:));

      flux_tot = sum(flux(:));

      fprintf(energyFile,'%e %e %e %e %e %e %e \n',t,dt,energy_tot,enstrophy_tot, ... 
              ZF_energy_tot,DW_energy_tot,flux_tot);

    end
    function plotfunc()
      % Go back in real space omega in real space for plotting
      diary off; %flush diary
      diary on;
      phi_curr = ifftz(phi_hat);
      w_curr=ifftz(phi_hat.*(-ksquare));
      n_curr=ifftz(n_hat);
      w = real(ifft2(w_curr));
      n = real(ifft2(n_curr)); 
      w = w(:,:,1);
      n = n(:,:,1);

      phi_curr = phi_curr(:,:,1);
      w_curr = w_curr(:,:,1);
      
      enstrophy = 0.5*w_curr.*conj(w_curr)/(NX*NY)^2;
      energy = 0.5*real(conj(phi_curr).*phi_curr.*ksquare(:,:,1))/(NX*NY)^2;                  

      if(padding)
         enstrophy=circshift(enstrophy,[NX_real/2,NY_real/2]); 
         enstrophy=enstrophy(2:NX_real,2:NY_real);
         energy=circshift(energy,[NX_real/2,NY_real/2]); 
         energy=energy(2:NX_real,2:NY_real);
      else
         enstrophy=circshift(enstrophy,[NX_real/2,NY_real/2]); 
         enstrophy=enstrophy(2:NX_real,2:NY_real);
         energy=circshift(energy,[NX_real/2,NY_real/2]); 
         energy=energy(2:NX_real,2:NY_real);
      end

      wlog=max(log10(enstrophy),-10);
      energylog=max(log10(energy),-10);
      m_w = squeeze(max(abs(w(:))));
      m_n = squeeze(max(abs(n(:))));
      set(0,'CurrentFigure',fig1);
      clf(fig1,'reset')
      cf=subplot(2,2,1);
      imagesc(LXnum,LYnum,w', [-m_w m_w]); axis equal tight; colorbar
      set(fig1.CurrentAxes,'Ydir','Normal')
      set(fig1.CurrentAxes,'Xdir','Normal')
      colormap(cf,c_maprb)
      title(sprintf('vorticity t=%.02f',t));
      xlabel('x');
      ylabel('y');
      cf=subplot(2,2,2);
      imagesc(LXnum,LYnum,n',[-m_n m_n]); axis equal tight; colorbar
      colormap(cf,c_maprb)
      title(sprintf('density t=%.02f',t));
      set(fig1.CurrentAxes,'Ydir','Normal')
      set(fig1.CurrentAxes,'Xdir','Normal')
      xlabel('x');
      ylabel('y');
      cf=subplot(2,2,3);
      imagesc(kxnum,kynum, energylog'); axis equal tight; colorbar
      colormap(cf,c_map)
      title('log10(Energy power spectrum)');    
      set(fig1.CurrentAxes,'Ydir','Normal')
      xlabel('kx');
      ylabel('ky');
      cf=subplot(2,2,4);
      imagesc(kxnum,kynum,wlog'); axis equal tight; colorbar
      colormap(cf,c_map)
      set(fig1.CurrentAxes,'Ydir','Normal')
      xlabel('kx');
      ylabel('ky');
      title('log10(vorticity/Enstrophy power spectrum)');    
      if(save_plots)
        save_binary_matrix(sprintf('%sdata/phi_%d.bin',basename,enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,phi(:,:,1));
        save_binary_matrix(sprintf('%sdata/w_%d.bin',basename,enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,w(:,:,1));
        save_binary_matrix(sprintf('%sdata/n_%d.bin',basename,enk),((1:NX)-0.5)*dx,((1:NY)-0.5)*dy,n(:,:,1));
        saveas(fig1,sprintf('%splots/fig_%d.png',basename,enk));
      end
      drawnow
      enk=enk+1;
    end

    function dump(filename)
        fileID=fopen(filename,'w');
        fwrite(fileID,t,'double');
        fwrite(fileID,dt,'double');
        fwrite(fileID,i,'int');
        fwrite(fileID,enk,'int');
        fwrite(fileID,trk,'int');
        fwrite(fileID,rsk,'int');
        fwrite(fileID,rek,'int');
        fwrite(fileID,dt,'double');
        fwrite(fileID,real(phi_hat),'double');fwrite(fileID,imag(phi_hat),'double');
        fwrite(fileID,real(n_hat),'double');fwrite(fileID,imag(n_hat),'double');
        fclose(fileID);
    end
end
