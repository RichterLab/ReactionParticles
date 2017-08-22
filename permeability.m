%permeability.m

%same as gaussian_correlated_new.m except written as a function

%Lx = length of domain in x
%Ly = length of domain in y
%var = variance of the lognormal field
%k_range = range of integers for wavelengths [-k_range*L,k_range*L]
%filenumber = number of file to save K field 
        %(i.e. if filenumber=1, we save perm1.mat)

function K = permeability(Lx,Ly,var,k_range,filenumber)
%K is lognormal

rng('shuffle')

%generate random field 
% Lx=10;
% Ly=1;

numgridy=5*Ly;
numgridx=5*Lx;

% Specify x-y grid
[x,y] = meshgrid(linspace(0,Lx,numgridx),linspace(0,Ly,numgridy));

y=flipud(y);

A=size(x);

n=A(1)*A(2);  %number of grid points %A(1)= number of rows, A(2) = number of columns

%var=1;  %variance of lognormal field

N=100;   %number of terms to sum over in generation of random field

dx=Lx/numgridx;
dy=Ly/numgridy;

%k_range=5; %integer range for possible wavelengths 
kx=randi([-k_range*Lx,k_range*Lx],1,N); %wavelength
ky=randi([-k_range*Ly,k_range*Ly],1,N);

phi=rand(1,N)*2*pi;  %random shift

sumK=zeros(A(1),A(2));       %normally disributed random field with just two components

for i=1:A(1) %sum over columns
    for j=1:A(2) %sum over rows
        for k=1:N   %N is the number of terms to sum over

            sumK(i,j)=sumK(i,j)+cos(2*pi/Lx*kx(k)*x(i,j)+2*pi/Ly*ky(k)*y(i,j)+phi(k));
            
        end
    end
end

K=exp(sqrt(2*var/N)*sumK); %lognormally distributed conductivity (2 components)
% K(:,N)=K(:,1);
% K(N,:)=K(1,:);

logK=log(K);

figure(1)
pcolor(x,y,logK)  %log(K)'
view([0 0 90])    %orientation of axes
axis([0 Lx 0 Ly])
title('logK')
colorbar
colormap jet
shading interp
% K=Kp;
daspect([1 1 1])  %data aspect ratio

% K_stat = zeros(1,A(1)*A(2));
% count=1;
% for i=1:A(1)
%     for j=1:A(2)
%         K_stat(count)=K(i,j);
%         count=count+1;
%     end
% end

% kk=1;
% average(kk)=mean(K_stat);
% variance(kk)=var(K_stat);

save(['perm',int2str(filenumber),'.mat'],'K')

end