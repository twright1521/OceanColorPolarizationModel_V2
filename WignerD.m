function [ WigDFunction ] = WignerD( m, n, orders, ScattAngle_deg )
%WignerD Computes the Wigner-D functions for the given m and n, over orders
%0:s, and all thetas. 
%
%   Inputs:
%    m,n    - the type of Wigner d-function required.   m=0,n=0 are equal
%             to the standard legendre polynomials.  
%    orders - The orders to return. (s) Should be smin:smax, where smin = 0
%             for reconstruction of F11 or F44 only.  smin = 2 for
%             reconstruction of F11+F44, F11-F44, F21, or F34
%  ScattAng - [degrees].  The angles over which you want the reconstruction
%             to be accurate.  
%
%  Output:
%  WigDFun  - Matrix where the rows pertain to the orders specified, and
%             the columns pertain to the angles. (See below).
%
%    Per Mischenko 2002, Appendix B
%
%         s
%        d  (theta)
%         m,n
%
% Output:
% WigDFunction(m,n,[s], [ScattAngle_deg]) =
%
%            . . . . .  cosd(ScattAngle_deg) . . . . . . . 
% row#      ________________________________________________
%   1      |  s = 0
%   2      |  s = 1
%   3      |  s = 2
%  ...     |   ...
% max(s)+1 |  s = max(s)
%
% Robert Foster

%Define x as the cosine of the angle
x = cosd(ScattAngle_deg);

%Get statistics about the input
smin = max(abs(m), abs(n));
smax = max(orders);

%Define the Psi function
if(n>=m)
    Psi = 1;
else
    Psi = (-1)^(m-n);
end



% WigDFunction(m,n,[s], [ScattAngle_deg]) =
%
%            . . . . .  cosd(ScattAngle_deg) . . . . . . . 
% row#      ________________________________________________
%   1      |  s = -1 (only for initializion)
%   2      |  s = 0
%   3      |  s = min(s)+1
%  ...     |   ...
% max(s)+1 |  s = max(s)
%

%Create the array to hold the output
%If max(s) = A, then we need A+2 rows, to include the 2 initial seeds
WigDFunction = nan(smax+2, length(x));

s2row = @(s)s+2;

%Define the initial value for smin-1 to be zero
WigDFunction(s2row(smin-1), :) = 0;

%Define the initial value for smin (s = 0 usually)
c1 = ...                               
    Psi * 2^(-smin) * sqrt(            factorial(2*smin) ...
                             / (factorial(abs(m-n))*factorial(abs(m+n))) );
c2 = (1-x).^(0.5*abs(m-n));
c3 = (1+x).^(0.5*abs(m+n));
WigDFunction(s2row(smin), :) = c1.*c2.*c3;




%Compute the Wigner d-function
%I don't even care to try to vectorize this anymore.
for t = 1:length(ScattAngle_deg)
   for s = smin:smax-1
       WigDFunction(s+3,t) = ...
          (    (2*s+1)*(s*(s+1)*x(t) - m*n)      * WigDFunction(s+2,t)  ...
             - (s+1)*sqrt(s^2-m^2)*sqrt(s^2-n^2) * WigDFunction(s+1,t)  ...
          ) ...
              ./( s*sqrt((s+1)^2 - m^2)*sqrt((s+1)^2 - n^2));

    %MATLAB IS NOT SMART ENOUGH TO REALIZE WHEN L'HOSPITAL'S RULE NEEDS TO
    %BE APPLIED.  [ (0x + 0)/0 can be = x ]
    % SO MANUALLY APPLY IT.  TOOK SO DAMN LONG TO FIGURE THIS OUT.
        if(isnan(WigDFunction(s+3,t)))
            WigDFunction(s+3,t) = x(t);
        end

          
   end
    
end

%Remove orders less than smin
WigDFunction(1:s2row(smin-1),:) = [];

end

