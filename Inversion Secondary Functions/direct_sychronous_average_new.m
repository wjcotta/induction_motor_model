function [out single] = direct_sychronous_average_new(data, wr, fs)

    L = length(data);
    md = mean(data);
            
    rps = wr/(2*pi);
    Period = 1/rps;
    
    new_fs = ceil(Period*fs)/Period;
    spp = round(Period*new_fs);
    new_L = round(L/fs*new_fs);
    
    data = interpft(data-md,new_L)+md;
    
    [avg, single] = sychronous_average(data);
     
    out = interpft(avg-md,L)+md;
    
    
    function [out h] = sychronous_average(meas)
        C = floor(length(meas)/spp);
        h = reshape(meas(1:C*spp),spp,C);
        h = mean(h,2);
        out = repmat(h,C,1);
        Lm = length(meas);
        out = [out; out(1:(Lm-(C*spp)))];
        
    end

end
