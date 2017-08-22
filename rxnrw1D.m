%1D reactive random walk model with particle killing method

function rxnrw1D(u,filename)

N = 10000;    %Initial Number of A (and B) Particles (Total number of particles = 2N)

D = 1e-5;     %Diffusion coefficient

L = 1;        %Length of the Domain
C0 = 1;       %Initial concentration
mp = C0*L/N;  %mass of each particle

k = 10;        %reaction rate

dt = 0.01;    %initial time step
Nsteps = 1e5; %number of timesteps
time = 0:dt:Nsteps*dt;

%rng('shuffle')
xA = L*rand(1,N);  %initial disribution of A particles at locations (xA,yA)
xB = L*rand(1,N);  %initial distribution of B particles

szxA = size(xA);

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

conc    = zeros(1,Nsteps+1);    %define concentration vector (add one to include t=0)
conc(1) = length(xA)*mp/L;           %initial concentration at t=0

meanU2 = zeros(1,Nsteps);
meanCA = zeros(1,Nsteps);

for kk = 1:Nsteps
    kk
    
    Pr = k*mp*dt;                 %probability of reaction given collocation

    %distance over which we consider that particles will react with one
    %another
    r = a*sqrt(2*D*dt);
    
    %reset particle velocities
    uA = zeros(szxA);
    uB = zeros(szxA);
    
    %Reaction step
    xdead = zeros(szxA);
    
    countA_grid = zeros(1,Ngridbox);
    countB_grid = zeros(1,Ngridbox);
    xBtemp = xB; %create a duplicate of xB for grid calculations
                 %(need this because xB may react before we get to it in
                 %the loop, but it is still needed for the <U^2>
                 %calculations)
    
    for ii = 1:length(xA)
        
        xAcurrent = xA(ii);
        
        %%%GRID FOR MEASURING CONCENCENTRATION QUANTITIES%%%%%%%%%%%%%%%%%
        rnddownxA = floor(xA(ii)/dx);
        rnddownxB = floor(xBtemp(ii)/dx);    
        
        %define indices of grid box in which the particle is located
        idx1xA = rnddownxA+1; %lower x index   %add 1 because index starts at 1, not 0
        idx1xB = rnddownxB+1; %lower x index
 
        countA_grid(idx1xA) = countA_grid(idx1xA)+1;
        countB_grid(idx1xB) = countB_grid(idx1xB)+1;
        
        %%%GRID FOR VELOCITY FIELD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rnddownxA_vel = floor(xA(ii)/dx_vel);
        rnddownxB_vel = floor(xBtemp(ii)/dx_vel);    
        
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
        
        uA(ii) = ((xgrid_vel(idx2xA_vel)-xAcurrent)*u(idx1xA_vel)+(xAcurrent-xgrid_vel(idx1xA_vel))*u(idx2xA_vel))...
                /(xgrid_vel(idx2xA_vel)-xgrid_vel(idx1xA_vel));

        uB(ii) = ((xgrid_vel(idx2xB_vel)-xB(ii))*u(idx1xB_vel)+(xB(ii)-xgrid_vel(idx1xB_vel))*u(idx2xB_vel))...
                /(xgrid_vel(idx2xB_vel)-xgrid_vel(idx1xB_vel));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        s1 = abs(xAcurrent-xB);
        s2 = abs(xAcurrent-xB-L); %account for periodicity, i.e. particles near an edge of
        s3 = abs(xAcurrent-xB+L); %the domain may react with particles near the opposite edge
        s = [s1; s2; s3];
        s = min(s); %only keep the shortest distance between each AB pair
        
        %calculate the probability of reaction of each AB particle pair given s
        Pf = Pr*1/sqrt(4*pi*(2*D)*dt)*exp(-s.^2/(4*(2*D)*dt));
        
        %generate a random number to determine if reaction will occur
        RP = Pf-rand(size(Pf));
        
        %identify the most probable of all possible reactions
        idxreact = find(RP==max(RP));
        
        if RP(idxreact)>0 %the AB particle pair reacts if RP > 0
            
            %AB react, so must be removed from the system according to
            %the reaction A+B->0
            xdead(ii) = 1; %indicate that this A particle reacted and needs to be removed
            xB(idxreact) = NaN; %use NaN as placeholder to later remove B particles that reacted
            
        end
        
    end   
    
    notdeadA = find(xdead==0); %indices of A particles that haven't reacted
    notdeadB = find(~isnan(xB)); %indices of B particles that haven't reacted
    xA = xA(notdeadA);
    xB = xB(notdeadB);
    uA = uA(notdeadA);
    uB = uB(notdeadB);
    
    szxA = size(xA);
    
    CA_grid = countA_grid*mp/dx;
    CB_grid = countB_grid*mp/dx;
    
    U = CA_grid-CB_grid;
    U2 = U.^2;
    
    meanCA(kk) = mean(CA_grid);
    meanU2(kk) = mean(U2);
    conc(1+kk) = length(xA)*mp/L;
    
    xA = xA+uA*dt+sqrt(2*D*dt)*randn(szxA);
    xB = xB+uB*dt+sqrt(2*D*dt)*randn(szxA);
    
    xA = mod(xA,L);  %periodic boundary condition
    xB = mod(xB,L);
    
    %stop simulation if all particles have reacted
    if length(xA)==0
        break
    end
    
end

figure(1)
hold on
plot(time,conc,'c.')
plot(time(1:end-1),meanCA,'k') %make sure these two are the exact same
set(gca,'Xscale','log','Yscale','log')
plot(time,C0./(1+k*C0*time),'r')

figure(2)
hold on
plot(time(1:end-1),meanU2,'k')
plot(time,conc,'c')
set(gca,'Xscale','log','Yscale','log')
legend('<U^2>','<CA>')

figure(3)
plot(log(time(30000:100000)),log(conc(30000:100000)))

save([filename,'.mat'],'conc','meanU2','meanCA','time')

end