function [ mie3InputFileName ] = writeMie3InputFile( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Define Default Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  D E F I N E  D E F A U L T  I N P U T  P A R A M E T E R S %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of different particles.  Must be between 1 and 4.
numParticles = 1;

%Expansion in generalized spherical functions (legendre) or Scat angle &
%F11, F21, F33, F34.   1 = Fourier-Legendre expansion, 0 = Scattering Angle
GSF = 1;

%Truncation of the Mie sum
mieSumTrunc = 1.0000000e-12;%originally 1.0000000e-12

%Cutoff value of the size distribution
distCutoff = 1.0000000e-11;

%Minimum scattering angle in degrees
minScattAngle_deg = 0;

%Maximum scattering angle in degrees
maxScattAngle_deg = 180;

%Scattering Angle step in degrees
scattAngleStep_deg = 1;

%THESE ARRAYS MUST BE OF LENGTH 4 EVEN IF USING <4 PARTICLES.  
%Wavelength (in micrometers)
lambda = [0.44, 2.25, 0.55, 0.45];

%Real part of the refractive index
Rem = [1.12, 1.42, 1.53, 1.33];

%Imaginary part of the refractive index (abs(Im(m)) ?) 
Imm = [0, 0, 0.0006, 0];

%   SizeDistIndex               PAR1        PAR2        PAR3                        
%  ---------------------------------------------------------                      
%   0 (no distribution)         r           -           -                           
%   1 (two parameter gamma)     alpha       b           -                           
%   2 (two parameter gamma)     reff        veff        -                           
%   3 (bimodal gamma)           reff1       reff2       veff                        
%   4 (log normal)              rg          sigma       -                           
%   5 (log normal)              reff        veff        -                           
%   6 (power law)               alpha       rmin        rmax                        
%   7 (modified gamma)          alpha       rc          gamma                       
%   8 (modified gamma)          alpha       b           gamma     

%Size distribution index
sizeDistIndex = [6, 5, 6, 8];

%Size Distribution Parameter 1 
%Applies to all size distribution indices
PAR1 = [4, 0.15, 4.005, 0.1];

%Size Distribution Parameter 2 
%Applies to all except 0 (no distribution)
%PUT 0 IF THIS PARAMETER DOESNT APPLY
PAR2 = [0.1, 0.2, 2.99, 8.9443];

%Size Distribution Parameter 3 
%(Applies to 3,6,7 & 8 (bimodal, power law and modified gamma)
%PUT 0 IF THIS PARAMETER DOESNT APPLY
PAR3 = [50, 0, 3, 0.5];

%Number of subintervals for r
rSubintervals = [50, 50, 40, 50];

%Number of Gauss points in sub
gaussPoints = [100, 100, 100, 100];

%% Process passed parameters and validate input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  V A L I D A T E  I N P U T  P A R A M E T E R S %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Instantiate the input parser
p = inputParser;

%Input Validation Functions
posIntLT4       = @(x) isscalar(x) && (x > 0) && (mod(x, 1)==0) && (x <=4);
posIntLT8       = @(x) all((x >= 0) & (mod(x, 1)==0) & (x <= 8));
posInt          = @(x) all((x >= 0) & (mod(x, 1)==0));
Eq01            = @(x) (x == 0) || (x == 1);
PosNumeric      = @(x) all(isnumeric(x) & (x >= 0));
Numeric         = @(x) all(isnumeric(x));
PosNumericLT180 = @(x) all(isnumeric(x) & (x >= 0) & (x >= 180));
Eq0_180         = @(x) (x >= 0) & (x <= 180);

%Add parameters to input parser
addParameter(p, 'numParticles', numParticles, posIntLT4);
addParameter(p, 'GSF', GSF, Eq01);
addParameter(p, 'distCutoff', distCutoff, PosNumeric);
addParameter(p, 'mieSumTrunc', mieSumTrunc, PosNumeric);
addParameter(p, 'minScattAngle_deg', minScattAngle_deg, Eq0_180);
addParameter(p, 'maxScattAngle_deg', maxScattAngle_deg, Eq0_180);
addParameter(p, 'scattAngleStep_deg', scattAngleStep_deg, PosNumericLT180);
addParameter(p, 'lambda', lambda, PosNumeric);
addParameter(p, 'Rem', Rem, PosNumeric);
addParameter(p, 'Imm', Imm, PosNumeric);
addParameter(p, 'sizeDistIndex', sizeDistIndex, posIntLT8);
addParameter(p, 'PAR1', PAR1, Numeric);
addParameter(p, 'PAR2', PAR2, Numeric);
addParameter(p, 'PAR3', PAR3, Numeric);
addParameter(p, 'rSubintervals', rSubintervals, posInt);
addParameter(p, 'gaussPoints', gaussPoints, posInt);

%Parse the results
parse(p, varargin{:});

%% Write the Mie input file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  W R I T E  M I E  I N P U T  F I L E %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Open the file
mie3InputFileName = 'mieCalcInput.txt';%_', datestr(now, 'HH_MM_SS_FFF'), '
FID = fopen(mie3InputFileName, 'wt');

%Write the file header
headerText = 'INPUT FOR THE MEERHOFF MIE PROGRAM VERSION 2.1                                 ';
fprintf(FID, ' %s\n\n', headerText);

%Write the number of particles
fprintf(FID, ' number of different particles (max = 4) : % 3u\n', p.Results.numParticles);

%Write the use GSF expansion bit
fprintf(FID, ' expansion in GSF (yes = 1, no = 0)      : % 3u\n', p.Results.GSF);

%Write section Divisions
fprintf(FID, ' ----------------------------------------- ---\n');
fprintf(FID, ' ----------------------------------------- ---.-------e---\n');

%Write Mie Sum Truncation
fprintf(FID, ' truncation of Mie sum                   : % 15.7e\n', p.Results.mieSumTrunc);

%Write Size distribution cutoff
fprintf(FID, ' cutoff value of the size distribution   : % 15.7e\n', p.Results.distCutoff);

%Write minimum scattering angle
fprintf(FID, ' minimum scattering angle in degrees     : % 15.7e\n', p.Results.minScattAngle_deg);

%Write maximum scattering angle
fprintf(FID, ' maximum scattering angle in degrees     : % 15.7e\n', p.Results.maxScattAngle_deg);

%Write scattering angle step
fprintf(FID, ' step in scattering angle in degrees     : % 15.7e\n', p.Results.scattAngleStep_deg);

%Write the next section text
divText =  ...
{' ----------------------------------------- ---.-------e---                      ';
' Definitions of size distributions are given below.                             ';
[];
'                            PAR1        PAR2        PAR3                        ';
' ---------------------------------------------------------                      ';
'  0 no distribution         r           -           -                           ';
'  1 two parameter gamma     alpha       b           -                           ';
'  2 two parameter gamma     reff        veff        -                           ';
'  3 bimodal gamma           reff1       reff2       veff                        ';
'  4 log normal              rg          sigma       -                           ';
'  5 log normal              reff        veff        -                           ';
'  6 power law               alpha       rmin        rmax                        ';
'  7 modified gamma          alpha       rc          gamma                       ';
'  8 modified gamma          alpha       b           gamma                       ';
[];
' PARTICLE NR. :                    1           2           3           4';
' ------------------------------ ---.------- ---.------- ---.------- ---.-------';};
for i = 1:length(divText)
    fprintf(FID, '%s\n', divText{i});
end

%Write the wavelength
fprintf(FID, ' %-31s', 'wavelength');
fprintf(FID, '%11.7f ', p.Results.lambda);

%Write the real part of index of refraction
fprintf(FID, '\n %-31s', 'Re(m)');
fprintf(FID, '%11.7f ', p.Results.Rem);

%Write imaginary part of index of refraction
fprintf(FID, '\n %-31s', 'abs(Im(m))');
fprintf(FID, '%11.7f ', p.Results.Imm);

%Write Size Distribution Index
fprintf(FID, '\n %-31s', 'index size distribution');
fprintf(FID, '% 3u         ', p.Results.sizeDistIndex);

%Write Number of subintervals for r
fprintf(FID, '\n %-31s', 'number of subintervals for r');
fprintf(FID, '% 3u         ', p.Results.rSubintervals);

%Write Number of Gauss points in sub
fprintf(FID, '\n %-31s', 'number of Gauss points in sub');
fprintf(FID, '% 3u         ', p.Results.gaussPoints);

%Write Parameter 1
fprintf(FID, '\n %-31s', 'PAR1');
fprintf(FID, '%11.7f ', p.Results.PAR1);

%Write Parameter 2
fprintf(FID, '\n %-31s', 'PAR2');
fprintf(FID, '%11.7f ', p.Results.PAR2);

%Write Parameter 3
fprintf(FID, '\n %-31s', 'PAR3');
fprintf(FID, '%11.7f ', p.Results.PAR3);

%Write footer text
footerText = ...
{' ------------------------------ ---.------- ---.------- ---.------- ---.-------';
[];
' DEFINITIONS OF SIZE DISTRIBUTIONS                                              ';
[];
' 0 NO DISTRIBUTION                                                              ';
'     single particle with radius r                                              ';
[];
' 1 TWO PARAMETER GAMMA                                                          ';
'     n(r) = C * r**alpha * exp(-b*r)                                            ';
[];
'     where     C = b**(alpha+1) / GAMMAFUNCTION(alpha+1)                        ';
[];
'     alpha must be larger than -1                                               ';
'     if alpha is less than zero, n(r) is singular at r = 0                      ';
[];
' 2 TWO PARAMETER GAMMA                                                          ';
'     n(r) = C * r**(1/veff-3) * exp(-r/(veff*reff))                             ';
[];
'     where     C = 1 / ((reff*veff)**(1/veff-2) * GAMMAFUNCTION(1/veff-2))      ';
[];
'     reff is the effective radius, reff must be positive                        ';
'     veff is the effective variance, veff must be positive and less than 0.5    ';
'     if veff is larger than 1/3, n(r) is singular at r = 0                      ';
[];
' 3 BIMODAL GAMMA                                                                ';
'     n(r) = 0.5 * ( TPG1 + TPG2 )                                               ';
[];
'     where TPG1 is a two parameter gamma distribution with reff = reff1         ';
'     and TPG2 is a two parameter gamma distribution with reff = reff2.          ';
'     TPG1 and TPG2 have the same effective variance veff.                       ';
'     reff1 must be positive, reff2 must be positive                             ';
'     veff must be positive and less than 0.5                                    ';
'     if veff is larger than 1/3, n(r) is singular at r = 0                      ';
[];
' 4 LOG NORMAL                                                                   ';
'     n(r) = C * (1/r) *  exp( -0.5*((log(r)-log(rg))/log(sigma))**2 )           ';
[];
'     where     C = 1 / (sqrt(2*pi)*log(sigma))                                  ';
[];
'     rg must be positive                                                        ';
'     sigma must be positive                                                     ';
[];
' 5 LOG NORMAL                                                                   ';
'     n(r) = C * (1/r) *  exp( -0.5*((log(r)-log(rg))/log(sigma))**2 )           ';
[];
'     where     rg = reff/(1+veff)**2.5 ,                                        ';
[];
'               log(sigma) = sqrt(log(1+veff))                                   ';
[];
'     and       C = 1 / (sqrt(2*pi)*log(sigma))                                  ';
[];
'     reff is the effective radius, reff must be positive                        ';
'     veff is the effective variance, veff must be positive                      ';
[];
' 6 POWER LAW                                                                    ';
'     n(r) =  C * r**(-alpha)        if rmin < r < rmax                          ';
'             0                      otherwise                                   ';
[];
'     where                                                                      ';
'     C = (alpha-1)*rmax**(alpha-1)/((rmax/rmin)**(alpha-1)-1)   if alpha not -1 ';
'         1/log(rmax/rmin)                                       if alpha is -1  ';
[];
'     rmin must be positive, rmax must be larger than rmin                       ';
[];
' 7 MODIFIED GAMMA                                                               ';
'     n(r) = C * r**alpha * exp(-b*r**gamma)                                     ';
[];
'     where   C = (gamma*b**((alpha+1)/gamma)) / GAMMAFUNCTION((alpha+1)/gamma)  ';
[];
'     and     b = alpha / (gamma*rc**gamma)                                      ';
[];
'     rc is the mode radius and must be positive                                 ';
'     alpha must be larger than -1                                               ';
'     gamma must be positive                                                     ';
'     if alpha is negative, n(r) is singular at r = 0                            ';
[];
' 8 MODIFIED GAMMA                                                               ';
'     n(r) = C * r**alpha * exp(-b*r**gamma)                                     ';
[];
'     where   C = (gamma*b**((alpha+1)/gamma)) / GAMMAFUNCTION((alpha+1)/gamma)  ';
[];
'     alpha must be larger than -1                                               ';
'     b must be positive                                                         ';
'     gamma must be positive                                                     ';
'     if alpha is negative, n(r) is singular at r = 0                            ';
[];
'   SPECIAL CASES OF 8 MODIFIED GAMMA (see Deirmendjian 1969, page 78)           ';
[];
'                                   Haze M       Haze L     Cloud C1             ';
'                       alpha      1.0000000    2.0000000   6.0000000            ';
'                       b          8.9443000   15.1186000   1.5000000            ';
'                       gamma      0.5000000    0.5000000   1.0000000            ';}; 
fprintf(FID, '\n');
for i = 1:length(footerText)
    fprintf(FID, '%s\n', footerText{i});
end

%Close the file handle
fclose(FID);
end

