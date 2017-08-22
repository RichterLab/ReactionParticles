xdead=zeros(size(xA));   %tracker to kill of A particles that have reacted

    
for jj=1:length(xA) %loop over A particles to determine which (if any) B particles they react with
       
        X=xA(jj)*ones(1,length(xB)); %the position of the A particle we are looking at
        Y=yA(jj)*ones(1,length(xB));
        
        %calculate minimum distance between this A and all B particles accounting for periodic
        %domain.
        s = [sqrt((X-xB).^2+(Y-yB).^2); sqrt((X-xB-Lx).^2+(Y-yB-Ly).^2); ...
             sqrt((X-xB-Lx).^2+(Y-yB).^2); sqrt((X-xB+Lx).^2+(Y-yB).^2); ...
		     sqrt((X-xB).^2+(Y-yB-Ly).^2); sqrt((X-xB).^2+(Y-yB+Ly).^2); ...
		     sqrt((X-xB-Lx).^2+(Y-yB+Ly).^2); sqrt((X-xB+Lx).^2+(Y-yB+Ly).^2); ...
		     sqrt((X-xB+Lx).^2+(Y-yB-Ly).^2)];
        s = min(s);
        
        %calculate the probability of reaction of each particle pair given
        %the distance between them
        Pf=Pr*1/(4*pi*(2*D)*dt)*exp(-s.^2/(4*(2*D)*dt));
        

        %generate a random number the size of each particle pair
        RP=Pf-rand(size(Pf));

        %identify the most probable of all possible reactions
        idxreact=find(RP==max(RP)); %index from list of B particles within range of A particle 

        %kill the B particle that partook in the most probable reaction
        if RP(idxreact)>0
            xdead(jj)=1; %indicate that the A particle has reacted
                       
            xB(idxreact)=NaN; %put NaN as a placeholder for the B particle that has reacted
            yB(idxreact)=NaN; %(will remove the B particle later)
        end
        
        %If they don't react, nothing happens. The particles stay in the
        %system and move by random walk in the next time step.
    
end
    
    %kill all particles that took place in a reaction
    notdead=find(xdead==0); %indices of A particles that haven't reacted    
    
    xA=xA(notdead);
    yA=yA(notdead);
    xB=xB(~isnan(xB));
    yB=yB(~isnan(yB));