function Ind = findIndex(PDF1, Goal, GetFutil)         % GetFutil = 1 if futility
L               = length(PDF1);
increment       = round(0.1*L);
Pos             = round(L/2);
S               = sum(PDF1(1:Pos));
increasing      = true;
if S>Goal
    increasing  = false;
end
if abs(sum(PDF1(1:1))) >= Goal
    increment   = 0;
    Pos         = 0;
end
if sum(PDF1) <= Goal
    increment   = 0;
    if GetFutil
        Pos	= 0;        % futility threshold
    else
        Pos = L;        % efficacy threshold
    end
end
while increment>1
    if increasing
        Pos = Pos+increment;
    else
        Pos = Pos-increment;
    end
    S = sum(PDF1(1:Pos));
    if increasing
        if S > Goal
            increment = floor(increment/2);
            increasing = false;
        end
    else
        if S < Goal
            increment = floor(increment/2);
            increasing = true;
        end
    end
end
Ind = Pos;