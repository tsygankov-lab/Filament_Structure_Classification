function [coeff,bnd_rec] = fourier_shape(bnd,bval,M,toreconst)

    xbnd = bnd(:,2); ybnd = bnd(:,1); nbnd = length(xbnd);
    coeff = [];
    bnd_rec = [];
        
    if max(xbnd)-min(xbnd)==0 || max(ybnd)-min(ybnd)==0
        return;
    end

    A = zeros(nbnd,2*M+1);
    A(:,1) = 1;
    for m=2:M+1
        A(:,m)   = cos((m-1)*(1:nbnd)/nbnd);
        A(:,m+M) = sin((m-1)*(1:nbnd)/nbnd);
    end
    
    xfft = A\xbnd;
    
    coeff = sqrt(xfft(2:M+1).^2+xfft(M+2:end).^2);  
    coeff = coeff/sum(coeff);
    
    if ~isempty(bval)
        vfft = A\bval;
        coeff2 = sqrt(vfft(2:M+1).^2+vfft(M+2:end).^2);  
        coeff2 = coeff2/sum(coeff2);
        coeff = [coeff; coeff2];
    end
    

    if toreconst
        
        yfft = A\ybnd;

        xrec = zeros(size(xbnd));
        yrec = zeros(size(xbnd));
        for n = 1:nbnd

            xrec(n) = xfft(1) + cos((1:M)*n/nbnd)*xfft(2:M+1) + sin((1:M)*n/nbnd)*xfft(M+2:end);
            yrec(n) = yfft(1) + cos((1:M)*n/nbnd)*yfft(2:M+1) + sin((1:M)*n/nbnd)*yfft(M+2:end);

        end
        bnd_rec = [yrec,xrec]; 
        
    end
    
    
    %xshp = [xfft(1); sqrt(xfft(2:M+1).^2+xfft(M+2:end).^2)];
    %xshp = xshp/sum(xshp);
    %yshp = [yfft(1); sqrt(yfft(2:M+1).^2+yfft(M+2:end).^2)];
    %yshp = yshp/sum(yshp);
    
    

