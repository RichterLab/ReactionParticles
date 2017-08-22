%rxnmass.m

%This is the 1D particle preserving method. Instead of removing particles
%that react we reduce their mass.

%clear
%clc
%close all

function rxnmass1D(u,filename)

N = 5e3; %Number of Particles
D = 1e-5; %diffusion coefficient
L = 1; %size of domain
k = 10; %reaction rate

%create initial condtion for particle locations
xA = L*rand(1,N);
xB = L*rand(1,N);

%mass vector for particles
mpA = 1/N*ones(size(xA));
mpB = 1/N*ones(size(xB));

dt = 0.1;  %time step

Nsteps = 2e3;  %number of time steps

a = 4; %factor multiplying diffusion distance - distance bigger than this do not react

%%%%GRID FOR CONCENTRATION QUANTITIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ngridx = 500*L; %number of grid points
Ngridbox = Ngridx-1; %number of grid cells
xgrid = linspace(0,L,Ngridx); %1D grid
dx = xgrid(2)-xgrid(1); %grid spacing

%%%%GRID FOR VELOCITY INTERPOLATION TO PARTICLES%%%%%%%%%%%%%%%%%%
Ngridx_vel = length(u); %number of grid points
Ngridbox_vel = Ngridx_vel-1; %number of grid cells
xgrid_vel = linspace(0,L,Ngridx_vel); %1D grid
dx_vel = xgrid_vel(2)-xgrid_vel(1); %grid spacing

%store vectors to measure total mass in domain at each time step
MA = zeros(1,Nsteps);
MB = zeros(1,Nsteps);

meanU2 = zeros(1,Nsteps);

%loop over timesteps
for kk=1:Nsteps

    kk

    mA_cell = zeros(1,Ngridbox); %matrix to store amount of mass in each grid cell
    mB_cell = zeros(1,Ngridbox);

    %loop over A particles
    for jj=1:length(xA)
        
        %set the A particle location and mass that will be used in this
        %loop
        xAcurrent = xA(jj);
        mpAcurrent = mpA(jj);
        
        %%%GRID FOR MEASURING CONCENCENTRATION QUANTITIES%%%%%%%%%%%%%%%%%
        rnddownA = floor(xA(jj)/dx);
        rnddownB = floor(xB(jj)/dx);    

        mA_cell(rnddownA+1) = mA_cell(rnddownA+1)+mpA(jj);
        mB_cell(rnddownB+1) = mB_cell(rnddownB+1)+mpB(jj);
        
        %%%GRID FOR VELOCITY FIELD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rnddownxA_vel = floor(xA(jj)/dx_vel);
        rnddownxB_vel = floor(xB(jj)/dx_vel);    
        
        %define indices of grid box in which the particle is located
        idx1xA_vel = rnddownxA_vel+1;           %lower x index   %add 1 because index starts at 1, not 0
        idx2xA_vel = rnddownxA_vel+2;           %upper x index   %could also use rndupxA+1
        idx1xB_vel = rnddownxB_vel+1;           %lower x index
        idx2xB_vel = rnddownxB_vel+2;           %upper x index
        
        %account for periodicity
        if idx2xA_vel==Ngridx_vel+1; 
            idx2xA_vel=1;
        end
        if idx2xB_vel==Ngridx_vel+1;
            idx2xB_vel=1;
        end
        
        %use bilinear interpolation to determine the velocity of the
        %particles        
        
        uA(jj) = ((xgrid_vel(idx2xA_vel)-xAcurrent)*u(idx1xA_vel)+(xAcurrent-xgrid_vel(idx1xA_vel))*u(idx2xA_vel))...
                /(xgrid_vel(idx2xA_vel)-xgrid_vel(idx1xA_vel));

        uB(jj) = ((xgrid_vel(idx2xB_vel)-xB(jj))*u(idx1xB_vel)+(xB(jj)-xgrid_vel(idx1xB_vel))*u(idx2xB_vel))...
                /(xgrid_vel(idx2xB_vel)-xgrid_vel(idx1xB_vel));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %distance over which we consider that particles will react with one
        %another
        r = a*sqrt(2*D*dt);
        
        %s is a vector storing distance from A particle to all other B
        %particles
        s = abs(xAcurrent-xB);
       
        %determine which B particles are within the reaction zone
        Bidx = find(s<r);     

        %loop over each particle pair within this distance
        for ii=1:length(Bidx)
            
            %distance between the two particles
            ss = s(Bidx(ii));   

            %change in mass of each particle due to the reaction occuring
            dm = k*mpAcurrent*mpB(Bidx(ii))/sqrt(8*pi*D*dt)*exp(-ss^2/(8*D*dt));
            
            %update particle masses both the A and he B particles
            mpAcurrent = mpAcurrent-dm*dt;
            mpB(Bidx(ii)) = mpB(Bidx(ii))-dm*dt;          

        end     

        %Now do the same thing accounting for periodicity at the lower and
        %upper boundaries - exact same as above, but checking distances
        %across lower and upper boundaries
        s = abs(xAcurrent-xB-L);     

        Bidx = find(s<r);     

        for ii = 1:length(Bidx)

            ss = s(Bidx(ii));
            
            dm = k*mpAcurrent*mpB(Bidx(ii))/sqrt(8*pi*D*dt)*exp(-ss^2/(8*D*dt));

            mpAcurrent = mpAcurrent-dm*dt;
            mpB(Bidx(ii)) = mpB(Bidx(ii))-dm*dt;           

        end

        
        s = abs(xAcurrent-xB+L);

        Bidx = find(s<r);

        
        for ii=1:length(Bidx)
          
            ss = s(Bidx(ii));

            dm = k*mpAcurrent*mpB(Bidx(ii))/sqrt(8*pi*D*dt)*exp(-ss^2/(8*D*dt));     

            mpAcurrent = mpAcurrent-dm*dt;
            mpB(Bidx(ii)) = mpB(Bidx(ii))-dm*dt;

        end

        %update the mass (in the ma vector) of the A particle that we have just looped over
        mpA(jj) = mpAcurrent;

    end

    %update total mass of A and B in the domain at the end of step kk
    MA(kk) = sum(mpA);
    MB(kk) = sum(mpB);

    CA = mA_cell/dx;
    CB = mB_cell/dx;

    U = CA-CB;
    U2 = U.^2;

    meanU2(kk) = mean(U2);

    %update particle positions by standard random walk
    xA = xA+uA*dt+sqrt(2*D*dt)*randn(size(xA));
    xB = xB+uB*dt+sqrt(2*D*dt)*randn(size(xA));

    %account for periodicity
    xA = mod(xA,L);
    xB = mod(xB,L);

end

%time vector for plotting
time = dt:dt:Nsteps*dt; 

CA = MA/L;
CB = MB/L; 


% figure(1)
% hold on
% plot(time,CA,'c.')
% plot(time,CB,'k') %make sure these two are the exact same
% set(gca,'Xscale','log','Yscale','log')
% plot(time,1./(1+time),'r')
% 
% figure(2)
% hold on
% plot(time,meanU2,'k')
% plot(time,CA,'c')
% set(gca,'Xscale','log','Yscale','log')
% legend('<U^2>','<CA>')
%
% figure(3)
% plot(log(time(1300:2000)),log(CA(1300:2000)))

save([filename,'.mat'],'CA','CB','meanU2','time')

end