function [equal] = isequaltol(A, B, tol)
%ISEQUALTOL 1 if max(abs(A-B))<tol, 0 else
%   This function checks two objects for almost equality up to a given
% tolerance tol.
equal = isequal(A, B);
if ~equal
    if iscell(A)
        if ~iscell(B)
            equal = 0;
            return;
        end
        allEqual=1;
        for k=1:length(A)
            allEqual = allEqual * isequaltol(A{k}, B{k}, tol);
        end
        equal = allEqual;
    else
        %inTol = ismembertol(A, B, tol);  
        inTol = abs(A-B)<tol;
        equal = min(min(min(min(min(min(min(inTol)))))));
    end
end
end

