% clc

alpha = 0.05;
gamma = [0.2, 0.4,0.25];
step = 0.0001;
for i = 1:3
    i
    if i==1
        xx = 0:step:20;
        phi = chi2pdf(xx,2);
        cummPhi = cumtrapz(xx,phi);
        
        idxA = find(1-cummPhi<alpha,1);
        
        A = xx(idxA)
        
        idxC = find(cummPhi>gamma(i),1);
        C = xx(idxC)
    else
        oldPhi = zeros(size(xx));
        oldPhi(idxC:idxA) =phi(idxC:idxA);
%         oldPhi = phi(idxC:idxA);
        phi = chi2pdf(xx,i+1);
        
%         c = filter2(phi,oldPhi);
        c = conv(oldPhi,chi2pdf(xx,i+3));
%         c =conv(oldPhi,chi2pdf(xx,i+1));
        
        %% Correcao?
%         phi = phi(1:2:size(phi,2))*step;
%         phi = phi(1:size(xx,2))*step;
        c = c*step;
        
        cummPhi = cumtrapz(xx,c(1:size(xx,2)));
        
%         [minValue,closestIndex] = min(abs())
        
        idxA = find(1-cummPhi<alpha,1);
%         idxA = find(1-cummPhi<alpha,1);
        A = xx(idxA)
        idxC = find(cummPhi>gamma(i),1);
        C = xx(idxC)
        
        figure(i)
        subplot(311)
        plot(xx,oldPhi)
        subplot(312)
        plot(xx,phi)
        subplot(313)
        plot(xx,c(1:size(xx,2)))
        
%         break;
    end
end

