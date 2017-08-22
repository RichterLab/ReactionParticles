function [u,v,P] = velocity_function(Lx,Ly,mean_velocity,filenumber)
% Lx=10;
% Ly=1;

load(['perm',int2str(filenumber),'.mat'])

M=5*Lx;
N=5*Ly;

deltax=Lx/M;
deltay=Ly/N;

deltaPx=1;
deltaPy=0;

Pfixed=1;   %specify pressure value (Pfixed) at a specific point P(1,N-1)

 
Kvec=reshape(K',[M*N,1]); %reshape so Kvec has same indices as P vector 
                       %(i.e. [11 21 31...M1; 12 22 32...M2; ...; 1N 2N
                       %       3N...MN] becomes [11; 21; 31; ...M1; 12; ...; MN])
A=spalloc(M*N,M*N,5*M*N);  %Equation we are solving is A*Pvec=B (Pvec=A\B)
B=zeros(M*N,1);
 
for i=1:M*N
    if mod(i,M)==1  %left boundary - P11,P12,P13,...,P1N
        
        if i==1 %corresponds to equation for P(1,1)
            
            A(i,i+M-2)=Kvec(i+M-2)/(deltax^2);  %P(M-1,1)
            A(i,i+1)=Kvec(i)/(deltax^2);        %P(2,1)
            A(i,i)=-((Kvec(i)+Kvec(i+M-2))/(deltax^2))-((Kvec(i)+Kvec(M*(N-2)+1))/(deltay^2)); %P(1,1)
            A(i,i+M*(N-2))=Kvec(i+M*(N-2))/(deltay^2); %P(1,N-1)
            A(i,i+M)=Kvec(i)/(deltay^2);               %P(1,2)
            B(i)=-Kvec(i+M-2)*deltaPx/(deltax^2)-Kvec(i+M*(N-2))*deltaPy/(deltay^2);
 
        elseif i==M*N-M+1 %corresponds to equation for P(1,N)
            
            A(i,i+M-2)=Kvec(i+M-2)/(deltax^2);  %P(M-1,N)
            A(i,i+1)=Kvec(i)/(deltax^2);        %P(2,N)
            A(i,i)=-((Kvec(i)+Kvec(i+M-2))/(deltax^2))-((Kvec(i)+Kvec(i-M))/(deltay^2)); %P(1,N)
            A(i,i-M)=Kvec(i-M)/(deltay^2);      %P(1,N-1)
            A(i,M+1)=Kvec(i)/(deltay^2);        %P(1,2)
            B(i)=-Kvec(i+M-2)*deltaPx/(deltax^2)+Kvec(i)*deltaPy/(deltay^2);
 
        else %left boundary interior points
            
            A(i,i+M-2)=Kvec(i+M-2)/(deltax^2);  %P(M-1,j)
            A(i,i)=-(Kvec(i)+Kvec(i+M-2))/(deltax^2)-(Kvec(i)+Kvec(i-M))/(deltay^2); %P(1,j)
            A(i,i+1)=Kvec(i)/(deltax^2);        %P(2,j)
            A(i,i-M)=Kvec(i-M)/(deltay^2);      %P(1,j-1)
            A(i,i+M)=Kvec(i)/(deltay^2);        %P(1,j+1)
            B(i)=-Kvec(i+M-2)*deltaPx/(deltax^2);
            
        end
    elseif mod(i,M)==0  %right boundary - PM1,PM2,PM3,...,PMN
        
        A(i,i)=-1;      %equation: deltaPx=P(1,j)-P(M,j)
        A(i,i-M+1)=1;
        B(i)=deltaPx;
        
    elseif i-M<0 && i-M>-M+1    %bottom boundary (P21, P31,...,P(M-1)1)
                                %no corners-they are already accounted for 
        A(i,i-1)=Kvec(i-1)/(deltax^2);                %P(i-1,1)
        A(i,i+1)=Kvec(i)/(deltax^2);                  %P(i+1,1)
        A(i,i)=-(Kvec(i)+Kvec(i-1))/(deltax^2)-(Kvec(i)+Kvec(i+M*(N-2)))/(deltay^2); %P(i,1)
        A(i,i+M*(N-2))=Kvec(i+M*(N-2))/(deltay^2);    %P(i,N-1)
        A(i,i+M)=Kvec(i)/(deltay^2);                  %P(i,2)
        B(i)=-Kvec(i+M*(N-2))*deltaPy/(deltay^2);
        
    elseif i+M>M*N+1 && i+M<M*N+M  %top boundary (P2N, P3N,...,P(M-1)N) -no corners
        
        A(i,i-M*(N-1))=1;      %equation: deltaPy=P(i,1)-P(i,N)
        A(i,i)=-1;
        B(i)=deltaPy;
 
%         A(i,i-1)=Kvec(i-1)/(deltax^2);                %P(i-1,N)
%         A(i,i+1)=Kvec(i)/(deltax^2);                  %P(i+1,N)
%         A(i,i)=(Kvec(i)+Kvec(i-1))/(deltax^2)+(Kvec(i)+Kvec(i-M))/(deltay^2); %P(i,N)
%         A(i,i-M)=Kvec(i-M)/(deltay^2);                %P(i,N-1)
%         A(i,i-M*(N-2))=Kvec(i)/(deltay^2);            %P(i,2)
%         B(i)=-Kvec(i)*deltaPy/(deltay^2);
        
    else  %all interior points of A
   
    A(i,i)=-(Kvec(i-1)/(deltax^2))-(Kvec(i-M)/(deltay^2))...
              -(Kvec(i)/(deltax^2))-(Kvec(i)/(deltay^2));          %P(i,j)
    A(i,i-1)=Kvec(i-1)/(deltax^2);                                 %P(i-1,j)
    A(i,i-M)=Kvec(i-M)/(deltay^2);                                 %P(i,j-1)
    A(i,i+1)=Kvec(i)/(deltax^2);                                   %P(i+1,j)
    A(i,i+M)=Kvec(i)/(deltay^2);                                   %P(i,j+1)
    
    end
end
 
%specify one point -choose P(1,N-1)
A(M*(N-2)+1,:)=zeros(1,M*N); %need to clear what is already in the row
A(M*(N-2)+1,M*(N-2)+1)=1;
B(M*(N-2)+1)=0;
B(M*(N-2)+1)=Pfixed;
 
%Find P by inverting A in the equation A*Pvec=B
Pvec=A\B;              %has indices [11; 21; 31; ...M1; 12; 22; 32; ...MN]
P=reshape(Pvec,[M,N])'; %has indices [11 12 13 ... 1N; 21 22 23 ..2N; 31; ...; MN]
 
u = zeros(N,M);
v = zeros(N,M);
 
%u and v are the velocities organized by index, u and v will later be organized
%by position
 
for i = 1:M-1
u(:,i)=-K(:,i).*(P(:,i+1)-P(:,i))/deltax;
end
 
u(:,M)=-K(:,M).*(P(:,2)-deltaPx-P(:,M))/deltax;  %right boundary (M+1=2 since periodic)
 
 
for j = 1:N-1
v(j,:)=-K(j,:).*(P(j+1,:)-P(j,:))/deltay;
end
 
v(N,:)=-K(N,:).*(P(2,:)-deltaPy-P(N,:))/deltay; %top boundary  (N+1=2 since periodic)
 
u=flipud(u);    %need to do flipud and take the transpose in order to
v=flipud(v);    %have the correct indices corresponding to position
                        %(11 is origin in bottom left, MN is top right)
                 % indices: [1N 2N 3N ... MN; ....; 12 22 32 ... M2; 11,21 31 ... M1]
 
% save u.mat u
% save v.mat v
meanu=mean(mean(u));
% meanv=mean(mean(v));
 
u = u./meanu.*mean_velocity;
v = v./meanu.*mean_velocity;

P=flipud(P);
 
save(['vel',int2str(filenumber),'.mat'],'u','v','P','-ascii');
 
figure(2)
pcolor(P);   
shading interp;       
colorbar;                  
colormap jet;           
title('Pressure')
daspect([1 1 1])
 
figure(3)
pcolor(u);
shading interp;
colorbar;
colormap jet;
title('u')
daspect([1 1 1]);
 
figure(4)
pcolor(v);
shading interp;
colorbar;
colormap jet;
title('v')
daspect([1 1 1]);

figure(5)
subplot(3,1,1)
pcolor(log(K));   
shading interp;       
colorbar;                  
colormap jet;           
title('log(\kappa)')
daspect([1 1 1])
set(gca, 'XTick', [], 'YTick',[]);
%text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,1,2)
pcolor(u);   
shading interp;       
colorbar;                  
colormap jet;           
title('u')
daspect([1 1 1])
set(gca, 'XTick', [], 'YTick',[]);
%text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,1,3)
pcolor(v);   
shading interp;       
colorbar;                  
colormap jet;           
title('v')
daspect([1 1 1])
set(gca, 'XTick', [], 'YTick',[]);
%text(0.02,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top')

end