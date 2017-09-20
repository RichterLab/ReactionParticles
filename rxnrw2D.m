%rxnrw.m

% clear
% clc
% close all

function rxnrw2D(u,v,filename)

N = 10;  % Initial Number of A (and B) Particles (Total number of particles = 2N)

D = 1e-5;     %Diffusion coefficient of A

Lx = 10;       %Length of the Domain
Ly = 1;
C0 = 1;       %Initial concentration
mp = C0*(Lx*Ly)/N;  %mass of each particle

k = 10;      %reaction rate

dt = 0.01;    %initial time step
epsilon = 1+0.025;    %rate at which timestep is increased
dtmax = 0.01;%8*pi*DA/k^2/mp^2/100;    %Maximum allowable timestep

Pr = k*mp*dt;     %probability of reaction
Pr_max = k*mp/8/pi/D;  %maximum probability of reaction (when s=0)
                       %need to make sure this number isn't too high

rng(1234)

xA = Lx*rand(1,N);  %initial disribution of A particles at locations (xA,yA)
xB = Lx*rand(1,N);  %initial distribution of B particles
yA = Ly*rand(1,N);  %initial disribution of A particles
yB = Ly*rand(1,N);  %initial distribution of B particles

xA(1) = 6.9234;
xA(2) = 4.3310;
xA(3) = 4.2018;
xA(4) = 1.4481;
xA(5) = 5.1285;
xA(6) = 8.0339;
xA(7) = 9.3188;
xA(8) = 6.4324;
xA(9) = 5.5442;
xA(10) = 1.1989;

yA(1) = 0.1467;
yA(2) = 0.1462;
yA(3) = 0.9986;
yA(4) = 0.0473;
yA(5) = 0.6285;
yA(6) = 0.6246;
yA(7) = 0.2868;
yA(8) = 0.2002;
yA(9) = 0.4078;
yA(10) = 0.9763;

xB(1) = 6.683352434689608;
xB(2) = 8.307203074836540;
xB(3) = 1.146952769896913;
xB(4) = 2.803967363212206;
xB(5) = 7.463246470539689;
xB(6) = 4.570490854194263;
xB(7) = 6.478712945735252;
xB(8) = 8.659819362606868;
xB(9) = 1.663233636699968;
xB(10) = 7.272637561834018;

yB(1) = 0.193934956352039;
yB(2) = 0.410571125591421;
yB(3) = 0.040404841745219;
yB(4) = 0.903354372755528;
yB(5) = 0.729667599112814;
yB(6) = 0.245990568047725;
yB(7) = 0.835649541816172;
yB(8) = 0.341178356375825;
yB(9) = 0.744876002240049;
yB(10) = 0.955182389042570;

Nsteps = 1;   %number of timesteps

%%%%%%%%%%CONCENTRATION GRID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ngridx = 2*Lx;     %number of gridpoints in x
Ngridy = 2*Ly;     %number of gridpoints in y

Ngridbox_x = Ngridx-1; %there are Ngrid-1 grid cells in each direction (Ngrid points)
Ngridbox_y = Ngridy-1;

%create grid - make (0,0) position bottom left corner
xgrid = repmat(linspace(0,Lx,Ngridx),Ngridy,1); %create grid of x positions (Ngrid x Ngrid)
ygrid = flipud(repmat(linspace(0,Ly,Ngridy),Ngridx,1)'); %grid of y positions

dx = xgrid(1,2) - xgrid(1,1);
dy = ygrid(1,1) - ygrid(2,1);

%%%%%%VELOCITY GRID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_vel = size(u);
Ngridx_vel = size_vel(2);     %number of gridpoints in x
Ngridy_vel = size_vel(1);     %number of gridpoints in y

%create grid - make (0,0) position bottom left corner
xgrid_vel = repmat(linspace(0,Lx,Ngridx_vel),Ngridy_vel,1); %create grid of x positions (Ngrid x Ngrid)
ygrid_vel = flipud(repmat(linspace(0,Ly,Ngridy_vel),Ngridx_vel,1)'); %grid of y positions

dx_vel = xgrid_vel(1,2) - xgrid_vel(1,1);
dy_vel = ygrid_vel(1,1) - ygrid_vel(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conc        = zeros(1,Nsteps+1);    %define concentration vector (add one to include t=0)
conc(1)     = length(xA)*mp/(Lx*Ly);           %initial concentration at t=0
time(1)     = 0;                    %define time vector, initial t=0;

meanU2      = zeros(1,Nsteps);
meanCA      = zeros(1,Nsteps);

for kk=1:Nsteps
    kk

    dt = min(dt*epsilon,dtmax);   %define timestep allowing it to increase to a maximum
    Pr = k*mp*dt;                 %probability of reaction given collocation

    P = 0.000001; %probability of reaction
    %r = sqrt(-8*D*dt*log(8*pi*D*dt*P/Pr)); %correct radius
    %radius of circle where particles at distance r have probability P of reacting
    szxA = size(xA);

    uA = zeros(szxA);
    uB = zeros(szxA);
    vA = zeros(szxA);
    vB = zeros(szxA);
    part = zeros(szxA);

    xBtemp = xB;
    yBtemp = yB;

    countA = zeros(Ngridy,Ngridx);
    countB = zeros(Ngridy,Ngridx);

    %bilinear interpolation to grid
    for ii=1:length(xA)

        xAnow = xA(ii);
        yAnow = yA(ii);
        xBnow = xBtemp(ii);
        yBnow = yBtemp(ii);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%CONCENTRATION GRID%%%%%%%%%
        %find out which gridbox each particle is in
        [xAnow, yAnow]
        rnddownxA = floor(xAnow/dx);
        rnddownyA = floor(yAnow/dy);
        rnddownxB = floor(xBnow/dx);
        rnddownyB = floor(yBnow/dy);

        %define indices of grid box in which the particle is located
        idx1xA = rnddownxA+1;           %lower x index   %add 1 because index starts at 1, not 0
        idx1yA = Ngridy-rnddownyA;      %lower y index
        idx1xB = rnddownxB+1;           %lower x index
        idx1yB = Ngridy-rnddownyB;      %lower y index

        %calculates the number of A and B particles in each gridbox
        [idx1yA, idx1xA]
        countA(idx1yA,idx1xA) = countA(idx1yA,idx1xA)+1; %add one because 0 is the first index
        countA(idx1yA,idx1xA)
        countB(idx1yB,idx1xB) = countB(idx1yB,idx1xB)+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%VELOCITY GRID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find out which gridbox each particle is in
        rnddownxA_vel = floor(xAnow/dx_vel);
        rnddownyA_vel = floor(yAnow/dy_vel);
        rnddownxB_vel = floor(xBnow/dx_vel);
        rnddownyB_vel = floor(yBnow/dy_vel);

        %define indices of grid box in which the particle is located
        idx1xA_vel = rnddownxA_vel+1;           %lower x index   %add 1 because index starts at 1, not 0
        idx2xA_vel = rnddownxA_vel+2;           %upper x index   %could also use rndupxA+1
        idx1yA_vel = Ngridy_vel-rnddownyA_vel;      %lower y index
        idx2yA_vel = Ngridy_vel-rnddownyA_vel-1;    %upper y index

        idx1xB_vel = rnddownxB_vel+1;           %lower x index
        idx2xB_vel = rnddownxB_vel+2;           %upper x index
        idx1yB_vel = Ngridy_vel-rnddownyB_vel;      %lower y index
        idx2yB_vel = Ngridy_vel-rnddownyB_vel-1;    %upper y index

        %account for periodicity
        if idx2xA_vel==Ngridx_vel+1;
            idx2xA_vel=1;
        end
        if idx2xB_vel==Ngridx_vel+1;
            idx2xB_vel=1;
        end
        if idx2yA_vel==0;
            idx2yA_vel=Ngridy_vel;
        end
        if idx2yB_vel==0;
            idx2yB_vel=Ngridy_vel;
        end

        %use bilinear interpolation to determine the velocity of the
        %particles
        uA(ii) = ((xgrid_vel(1,idx2xA_vel)-xAnow)*(yAnow-ygrid_vel(idx1yA_vel,1))*u(idx2yA_vel,idx1xA_vel)...
                +(xgrid_vel(1,idx2xA_vel)-xAnow)*(ygrid_vel(idx2yA_vel,1)-yAnow)*u(idx1yA_vel,idx1xA_vel)...
                +(xAnow-xgrid_vel(1,idx1xA_vel))*(ygrid_vel(idx2yA_vel,1)-yAnow)*u(idx1yA_vel,idx2xA_vel)...
                +(xAnow-xgrid_vel(1,idx1xA_vel))*(yAnow-ygrid_vel(idx1yA_vel,1))*u(idx2yA_vel,idx2xA_vel))...
                /((xgrid_vel(1,idx2xA_vel)-xgrid_vel(1,idx1xA_vel))*(ygrid_vel(idx2yA_vel,1)-ygrid_vel(idx1yA_vel,1)));

        uB(ii) = ((xgrid_vel(1,idx2xB_vel)-xBnow)*(yBnow-ygrid_vel(idx1yB_vel,1))*u(idx2yB_vel,idx1xB_vel)...
                +(xgrid_vel(1,idx2xB_vel)-xBnow)*(ygrid_vel(idx2yB_vel,1)-yBnow)*u(idx1yB_vel,idx1xB_vel)...
                +(xBnow-xgrid_vel(1,idx1xB_vel))*(ygrid_vel(idx2yB_vel,1)-yBnow)*u(idx1yB_vel,idx2xB_vel)...
                +(xBnow-xgrid_vel(1,idx1xB_vel))*(yBnow-ygrid_vel(idx1yB_vel,1))*u(idx2yB_vel,idx2xB_vel))...
                /((xgrid_vel(1,idx2xB_vel)-xgrid_vel(1,idx1xB_vel))*(ygrid_vel(idx2yB_vel,1)-ygrid_vel(idx1yB_vel,1)));

        vA(ii) = ((xgrid_vel(1,idx2xA_vel)-xAnow)*(yAnow-ygrid_vel(idx1yA_vel,1))*v(idx2yA_vel,idx1xA_vel)...
                +(xgrid_vel(1,idx2xA_vel)-xAnow)*(ygrid_vel(idx2yA_vel,1)-yAnow)*v(idx1yA_vel,idx1xA_vel)...
                +(xAnow-xgrid_vel(1,idx1xA_vel))*(ygrid_vel(idx2yA_vel,1)-yAnow)*v(idx1yA_vel,idx2xA_vel)...
                +(xAnow-xgrid_vel(1,idx1xA_vel))*(yAnow-ygrid_vel(idx1yA_vel,1))*v(idx2yA_vel,idx2xA_vel))...
                /((xgrid_vel(1,idx2xA_vel)-xgrid_vel(1,idx1xA_vel))*(ygrid_vel(idx2yA_vel,1)-ygrid_vel(idx1yA_vel,1)));

        vB(ii) = ((xgrid_vel(1,idx2xB_vel)-xBnow)*(yBnow-ygrid_vel(idx1yB_vel,1))*v(idx2yB_vel,idx1xB_vel)...
                +(xgrid_vel(1,idx2xB_vel)-xBnow)*(ygrid_vel(idx2yB_vel,1)-yBnow)*v(idx1yB_vel,idx1xB_vel)...
                +(xBnow-xgrid_vel(1,idx1xB_vel))*(ygrid_vel(idx2yB_vel,1)-yBnow)*v(idx1yB_vel,idx2xB_vel)...
                +(xBnow-xgrid_vel(1,idx1xB_vel))*(yBnow-ygrid_vel(idx1yB_vel,1))*v(idx2yB_vel,idx2xB_vel))...
                /((xgrid_vel(1,idx2xB_vel)-xgrid_vel(1,idx1xB_vel))*(ygrid_vel(idx2yB_vel,1)-ygrid_vel(idx1yB_vel,1)));

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    countA
    countB

   
    %remove row and column of zeros in countA and countB
    countA = countA(2:end,1:end-1);
    countB = countB(2:end,1:end-1);
    
    
    CA = countA*mp/(dx*dy);
    CB = countB*mp/(dx*dy);

    U = CA-CB;
    U2 = U.^2;

    meanU2(kk) = mean(mean(U2));
    meanCA(kk) = mean(mean(CA));

    %Random Walk part - update the location of particles x and y by a Brownian motion
    part = randn(szxA);
    xA = xA+uA*dt+sqrt(2*D*dt)*part;

    part = randn(szxA);
    xB = xB+uB*dt+sqrt(2*D*dt)*part;

    part = randn(szxA);
    yA = yA+vA*dt+sqrt(2*D*dt)*part;

    part = randn(szxA);
    yB = yB+vB*dt+sqrt(2*D*dt)*part;

    xA = mod(xA,Lx);  %periodic boundary
    xB = mod(xB,Lx);
    yA = mod(yA,Ly);
    yB = mod(yB,Ly);

    posA = [xA;yA]; %gives position matrix of x in row 1, y in row 2
    posB = [xB;yB];

    xyA = transpose(posA); %first column x, second column y
    xyB = transpose(posB);

    %stop simulation if all particles have reacted
    if isempty(xyA)==1
        break
    end

    %RXN
    rxn_cycle2D    %uses kd-tree to search particles more quickly


    %calculate the concentration at a given time
    conc(kk+1) = length(xA)*mp/(Lx*Ly);  %conc = num part/domain volume = num part/(L^2) = num part
    time(kk+1) = time(kk)+dt;

end

%renormalize concentration to 1 at t=0;
conc=conc/N;     %normalize by initial conc

save([filename,'.mat'],'conc','meanU2','meanCA','time');

%plot concentration against time

 %figure(1)
 %plot(time,conc(1:length(time)),'k')
 %set(gca,'Xscale','log','Yscale','log')
 %hold on
 %plot(time,1./(1+k*C0*time),'r')
 %legend('rxnrw','well-mixed')
 %xlabel('time')
 %ylabel('concentration')
end