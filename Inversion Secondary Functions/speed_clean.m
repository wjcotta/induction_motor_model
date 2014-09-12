function Speed = speed_clean(counts, divisor, cpr)

range = 0.5; % in rad/sec

if nargin < 3
    cpr = 2500;
end

if size(counts,2) == 2
    
    A = counts(:,1);
    B = counts(:,2);
    L = length(A);
    SpeedA = 2*pi*(48000000*(divisor/2)./(cpr*(A)));
    SpeedB = 2*pi*(48000000*(divisor/2)./(cpr*(B)));
    mA = medfilt1(SpeedA,10,1e4); % medfilt1 is a one dimensional median 
    % filter
    mB = medfilt1(SpeedB,10,1e4);
    indA = find(abs(SpeedA-mA)>range); % find returns the indices of the 
    % specified elements
    indB = find(abs(SpeedB-mB)>range);
    
    for i = 1:length(indA)
        if ~ismember(indA(i),indB)
            A(indA(i)) = B(indA(i));
        elseif ~ismember(indA(i)+1,indB) && (indA(i)+1)<=L
            A(indA(i)) = B(indA(i)+1);
        elseif ~ismember(indA(i)-1,indB) && (indA(i)-1)>0
            A(indA(i)) = B(indA(i)-1);
        elseif ~ismember(indA(i)+2,indB) && (indA(i)+2)<=L
            A(indA(i)) = B(indA(i)+2);
        elseif ~ismember(indA(i)-2,indB) && (indA(i)-2)>0
            A(indA(i)) = B(indA(i)-2);
        elseif ~ismember(indA(i)+3,indB) && (indA(i)+3)<=L
            A(indA(i)) = B(indA(i)+3);
        elseif ~ismember(indA(i)-3,indB) && (indA(i)-3)>0
            A(indA(i)) = B(indA(i)-3);
        end
    end
    
    for i = 1:length(indB)
        if ~ismember(indB(i),indA)
            B(indB(i)) = A(indB(i));
        elseif ~ismember(indB(i)+1,indA) && (indB(i)+1)<=L
            B(indB(i)) = A(indB(i)+1);
        elseif ~ismember(indB(i)-1,indA) && (indB(i)-1)>0
            B(indB(i)) = A(indB(i)-1);
        elseif ~ismember(indB(i)+2,indA) && (indB(i)+2)<=L
            B(indB(i)) = A(indB(i)+2);
        elseif ~ismember(indB(i)-2,indA) && (indB(i)-2)>0
            B(indB(i)) = A(indB(i)-2);
        elseif ~ismember(indB(i)+3,indA) && (indB(i)+3)<=L
            B(indB(i)) = A(indB(i)+3);
        elseif ~ismember(indB(i)-3,indA) && (indB(i)-3)>0
            B(indB(i)) = A(indB(i)-3);
        end
    end
    
    clean = A+B;
    Speed = 2*pi*(48000000*divisor./(cpr*(clean)));
else % counts only one column.
    

    A = counts(:,1);
    
    Speed = 2*pi*(48000000*(divisor)./(cpr*(A)));
    mA = medfilt1(Speed,5,1e4);
    indA = find(abs(Speed-mA)>range);
    
    for i = 1:length(indA)
        Speed(indA(i)) = mA(indA(i));
    end
   
end




end
