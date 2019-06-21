function [  a_seawater, b_seawater, c_seawater,                 ...
            n_SeaWater, excessScattPercent, h2o_density_kgm3,   ...
            a_seawater_std, b_seawater_std, c_seawater_std]     ...
            = SeaWaterIOPs( Wavelength_nm, varargin )
%SeaWaterIOPs Computes the scattering, absorption, and attenuation of water
%             given its Salinity, Temperature, Depth, Depolarization ratio,
%             and wavelength.  Also returns the refractive index of the
%             water, the excess scattering due to salinity, the water
%             density, and uncertainty estimates for scattering, absorption
%             and attenuation spectra.
%
%   REQUIRED INPUTS:
%       Wavelength_nm           The wavelength in nanometers.  Can be an
%                               array of wavelengths.
%
%   OPTIONAL INPUTS:
%       Salinity_ppt            The Salinity of the desired water in parts
%                               per thousand.  If not given, assumes a
%                               value of 32 ppt.
%       Temperature_C           The temperature of the water, expressed in
%                               degrees Celsius.  If not given, assumes the
%                               temperature is 20C.
%       DepolRatio              The Rayleigh Depolarization ratio for the
%                               water.  Morel determined this value to be
%                               0.09, however more recent studies
%                               determined 0.039 is best.  Defaults to
%                               0.039 if not given.
%       Depth_m                 The depth in meters below the sea surface
%                               of the water in which you are interested.
%                               Affects the results mostly through the
%                               change in density.  Defaults to 0 meters if
%                               not given. Calculation at > 1 atm may be
%                               unreliable.
%       RefrIdx_Method          The method by which the refractive index of
%                               the water is calculated.  The options are
%                               "QuanFry", or "Eisenberg".  Quan & Fry is
%                               very accurate for shallow depths.  The
%                               default is "QuanFry", unless the Depth_m
%                               given is > 100m.
%       Density_Method          The method by which the density of water is
%                               calculated.  Options are "UNESCO",
%                               "BryanCox", or "Mamayev".  UNESCO is the
%                               default.  Mamayev is only valid for shallow
%                               depths.
%       BetaT_Method            The method by which the isothermal
%                               compressibility of water is calculated.
%                               Option are "Kell" and "Lepple".  Kell is
%                               the default. Both methods are only valid
%                               for 1 atm.  Results could be slightly
%                               inaccurate at large depths.
%
%   OUTPUTS:
%       a_seawater              The absorption coefficient of the described
%                               seawater, at the given wavelengths.  The
%                               units are [1/m].
%       b_seawater              The scattering coefficient of the described
%                               seawater, at the given wavelengths.  The
%                               units are [1/m].
%       c_seawater              The attenuation coefficient of the 
%                               described seawater, at the given 
%                               wavelengths.  Equal to the sum of
%                               a_seawater and b_seawater. The units are 
%                               [1/m].
%       n_seawater              The refractive index of the described sea
%                               water, at the given wavelengths, calculated
%                               using the "RefrIdx_Method" described above.
%                               The index returned is w.r.t. vacuum.
%       excessScattPercent      The percentage increase in scattering
%                               of the described water over the same water
%                               without any salt ions. (Salitinty_ppt = 0),
%                               excessScattPercent is given as a decimal
%                               value [0 <= excessScattPercent <= 1].
%       h20_density_kgm3        The density of the described seawater, in
%                               units of [kg/m^3].
%       a_seawater_std          A 1 sigma confidence uncertainty estimate
%                               for the absorption spectra at the given
%                               wavelengths, in units of [1/m].
%       b_seawater_std          A 1 sigma confidence uncertainty estimate
%                               for the scattering spectra at the given
%                               wavelengths, in units of [1/m].
%       c_seawater_std          A 1 sigma confidence uncertainty estimate
%                               for the attenuation spectra at the given
%                               wavelengths, in units of [1/m].
%                               
%                               
%
%   USAGE:
%
%   Minimal use: (assumes 32ppt salinity, 20C, 0.039 depolarization ratio)
%   [a, b, c, n, excessScatt, rho, a_std, b_std, c_std] =               ...
%                                       SeaWaterIOPs(400:20:700);
%
%   Typical Use: (define salinity, temperature, depolatization ratio)
%   [a, b, c, n, excessScatt, rho, a_std, b_std, c_std] =               ...
%                          SeaWaterIOPs(400:20:700,                     ...
%                                      'Salinity_ppt',    32,           ...
%                                      'Temperature_C',   20,           ...
%                                      'DepolRatio',      0.039);
%
%   With all options specified:
%   [a, b, c, n, excessScatt, rho, a_std, b_std, c_std] =               ...
%                          SeaWaterIOPs(400:20:700,                     ...
%                                      'Salinity_ppt',    32,           ...
%                                      'Temperature_C',   20,           ...
%                                      'DepolRatio',      0.039,        ...
%                                      'Depth_m',         0,            ...
%                                      'RefrIdx_Method',  'QuanFry',    ...
%                                      'BetaT_Method',    'Kell',       ...
%                                      'Density_Method',  'UNESCO');
%
%   REFERENCES:
%       Bryan, K. and M. D. Cox (1972). "An Approximate Equation of State
%       for Numerical Models of Ocean Circulation." Journal of Physical
%       Oceanography 2(4): 510-514.
% 
%       Buiteveld, H., et al. (1994). Optical properties of pure water.
% 
%       Cabannes, J. (1920). "Relation entre le degré de polarisation et
%       l'intensité de la lumière diffusée par des molécules anisotropes.
%       Nouvelle détermination de la constante d'Avogadro." J. phys. radium
%       1(5):129-142.
% 	
%       Ciddor, P. E. (1996). "Refractive index of air: new equations for
%       the visible and near infrared." Applied Optics 35(9): 1566-1573.
%
%       Eisenberg, H. (1965). "Equation for the Refractive Index of Water."
%       The Journal of Chemical Physics 43(11): 3887-3892.
% 	
%       Evtyushenkov, A. and Y. F. Kiyachenko (1982). "Determination of the
%       dependence of liquid refractive index on pressure and temperature."
%       Optics and Spectroscopy 52: 56-58.
% 	
%       Fofonoff, N. P. and R. C. Millard (1983). Algorithms for
%       computation of fundamental properties of seawater. Unesco technical
%       papers in marine science. 44.
% 
%       Jonasz, M. and G. R. Fournier (2007). Light Scattering By Particles
%       In Water. Amsterdam, Academic Press.
% 
%       Kell, G. S. (1970). "Isothermal compressibility of liquid water at
%       1 atm." Journal of Chemical and Engineering Data 15(1): 119-122.
% 	
%       Kratohvil, J. P., et al. (1965). "Light Scattering by Pure Water."
%       The Journal of Chemical Physics 43(3): 914-921.
% 	
%       Lepple, F. K. and F. J. Millero (1971). "The isothermal
%       compressibility of seawater near one atmosphere." Deep Sea Research
%       and Oceanographic Abstracts 18(12): 1233-1254.
% 
%       Mamayev, O. I. (1975). Temperature-Salinity Analysis of World Ocean
%       Waters. Amsterdam, Netherlands, Elsevier Science.
% 	
%       Millero, F. J., et al. (1981). Background papers and supporting
%       data on the International Equation of State of Seawater, 1980.
%       Unesco technical papers in marine science, Unesco/ICES/SCOR/IAPSO
%       Joint Panel on Oceanographic Tables and Standards: 192.
% 	
%       Millero, F. J. and W. H. Leung (1976). "Thermodynamics of seawater
%       at one atmosphere." Am. J. Sci. 276(9): 1035-1077.
% 	
%       Morel, A. (1974). "Optical properties of pure water and pure sea
%       water." Optical aspects of oceanography 1: 1-24.
% 
%       O'Connor, C. L. and J. P. Schlupf (1967). "Brillouin Scattering in
%       Water: The Landau—Placzek Ratio." The Journal of Chemical Physics
%       47(1): 31-38.
%
%       Pegau, W. S., et al. (1997). "Absorption and attenuation of visible
%       and near-infrared light in water: dependence on temperature and
%       salinity." Applied Optics 36(24): 6035-6046.
% 
%       Pope, R. M. and E. S. Fry (1997). "Absorption spectrum (380?700 nm)
%       of pure water. II. Integrating cavity measurements." Applied Optics
%       36(33): 8710-8723.
%
%       Proutiere, A., et al. (1992). "Refractive index and density
%       variations in pure liquids: a new theoretical relation." The
%       Journal of Physical Chemistry 96(8): 3485-3489.
%	
%       Rayleigh, L. (1920). "A re-examination of the light scattered by
%       gases in respect of polarisation. I. Experiments on the common
%       gases." Proceedings of the Royal Society of London. Series A,
%       Containing Papers of a Mathematical and Physical Character:
%       435-450.
% 	
%       Zhang, X. and L. Hu (2009). "Estimating scattering of pure water
%       from density fluctuation of the refractive index." Optics Express
%       17(3): 1671-1678.
%
%       Zhang, X., et al. (2009). "Scattering by pure seawater: Effect of
%       salinity." Optics Express 17(7): 5698-5710.
% 
%       Zhang, X., et al. (2009). "Scattering by solutions of major sea
%       salts." Optics Express 17(22): 19580-19585.
% 
% Written by:
% Robert Foster
% NOAA-CREST Optical Remote Sensing Lab
% Grove School of Engineering, Rm T-553
% The City College of New York
% Email: rfoster01@citymail.cuny.edu
%
% UPDATES
% Version 1.0:  5/19/2015
% Version 2.0:  5/28/2015
%        - Changed Riso and Rcf computation method to Zhang, 2009a,b,c
%        - Zhang method required Proutiere for rho*(dn^2/drho)
%        - Computed dn^2dS from Quan & Fry Formula
%        - Eisenberg Method seems not to work well for seawater.
% Version 3.0   6/4/2015
%        - Compute uncertainty estimates for aw, bw, cw.
%

%==========================================================================
%%%%%%%%%%%%%%%%%%% V A L I D A T E  I N P U T %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%Instantiate the input parser
p = inputParser;

%Input Validation Functions
wavelenLimits   = ...  (400-700)from Quan & Fry 1995 refr. index formula
    @(x) all(isnumeric(x) & (x >= 350) & (x <= 950) & isvector(x)); 
salinityLimits  = ...  35 is from Quan & Fry 1995
    @(x) isscalar(x) && isnumeric(x) && (x >= 0) && (x <= 40);  
tempCLimits     = ...  Limited by Bryan & Cox 1972 Density Formula
    @(x) isscalar(x) && isnumeric(x) &&(x >= -2) && (x <= 30);   
depthLimits     = ...  Limited by Bryan & Cox 1972 Density Formula
    @(x) isscalar(x) && isnumeric(x) &&(x >= 0) && (x <= 6000);
BetaTChoices    = ... Tried both, I like Kell better. Both good.
    @(y) strcmp('Kell',y) || strcmp('Lepple',y);
RefrIdxChoices  = ... Quan & Fry is better, but doesnt deal with density
    @(y) strcmpi('QuanFry',y) || strcmpi('Eisenberg',y);
DensityChoices  = ... Bryan&Cox 1972, Mamayev 1975 or UNESCO
    @(y) strcmpi('BryanCox',y) || strcmpi('Mamayev',y) ...
      || strcmpi('UNESCO',y);
DepolLimits     = ... Depolarization ratio should not be more than 0.2
    @(x) isnumeric(x) && (x >= 0) && (x <= 0.2);

%Parameters    Name                 DefaultValue    Validation Function
addRequired(p, 'Wavelength_nm',                     wavelenLimits   );
addOptional(p, 'Salinity_ppt',      32,             salinityLimits  );
addOptional(p, 'Temperature_C',     20,             tempCLimits     );
addOptional(p, 'Depth_m',           0,              depthLimits     );
addOptional(p, 'DepolRatio',        0.039,          DepolLimits     );
addOptional(p, 'BetaT_Method',      'Kell',         BetaTChoices    );
addOptional(p, 'RefrIdx_Method',    'QuanFry',      RefrIdxChoices  );
addOptional(p, 'Density_Method',    'UNESCO',       DensityChoices  );

%Parse the results
try    parse(p, Wavelength_nm, varargin{:});
catch  err
       fprintf(2, 'Input Failed Validation.\n');
       fprintf(2, '%s\n', err.identifier);
       rethrow(err);
end

%Assign validated inputs
wavelength_nm   = p.Results.Wavelength_nm;
salinity_ppt    = p.Results.Salinity_ppt;
temperature_C   = p.Results.Temperature_C;
depth_m         = p.Results.Depth_m;
delta           = p.Results.DepolRatio;
BetaT_Method    = p.Results.BetaT_Method;
RefrIdx_Method  = p.Results.RefrIdx_Method;
Density_Method  = p.Results.Density_Method;

%NOTE: Quan & Fry's formula is only valid for 1 atm. Eisenbergs formula is
%dubious though, so use QF for depths of 100m or less (guessing)
if(depth_m > 100 && ~strcmpi(RefrIdx_Method, 'Eisenberg'))
    RefrIdx_Method = 'Eisenberg';
    warning(['Quan & Frys Refractive Index formula is only valid ', ...
        'for ~1 atm of pressure.  Using Eisenbergs Method.']);
end

%NOTE: Mamayev's Density formula is only valid for 1 atm. If the depth >100
%meters, use UNESCO method.
if(depth_m > 100 && strcmpi(Density_Method, 'Mamayev'))
    Density_Method = 'UNESCO';
    warning(['Mamayevs density formula is only valid ', ...
        'for ~1 atm of pressure.  Using UNESCO Method.']);
end

%Warn if the wavelength is outside of the range for Quan and Fry
% if(strcmpi(RefrIdx_Method, 'QuanFry') && ...
%                 (any(wavelength_nm < 400) || any(wavelength_nm > 700)))
%     warning(['Quan & Fry Refractive Index Method may not be accurate ',...
%         'less than 400nm or greater than 700nm.']);
% end

%Warn if the salinity is outside of the range for Quan and Fry
if(strcmpi(RefrIdx_Method, 'QuanFry') && (salinity_ppt > 35))
    warning(['Quan & Fry Refractive Index Method may not be accurate ',...
        'for Salinity greater than 35ppt.']);
end

%If the wavelength input array is a column
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%% C O N V E R T  U N I T S %%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================

%Convert the wavelength to meters and micrometers
wavelength_m = wavelength_nm .* 10^-9;
wavelength_um = wavelength_nm ./ 1000;

%Convert the temperature to Kelvin
temperature_K = temperature_C + 273.15;

%Convert salinity from parts per thousand to grams/gram
S_gpg = salinity_ppt / 1000;     % [g/g]

%==========================================================================
%%%%%%%%%%%%%%%%%%% D E F I N E  C O N S T A N T S %%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%Boltzmann Constant
%Units = J/K = (kg*m^2)/(K*s^2) =  Pa*(m^3/K)
Kb = 1.3806488 * 10^-23;  

% Avogadro's number [1/mol]
Na = 6.0221412927 * 10^23; 

%Molecular weight of water  [g/mol]
MolWeight_h20 = 18.01528;

%==========================================================================
%%%% C O M P U T E  I S O T H E R M A L  C O M P R E S S I B I L I T Y %%%%
%==========================================================================

%      ---------------------[ Method #1 ]-----------------------
if(strcmpi(BetaT_Method, 'Lepple'))
    %From Lepple and Millero 1971  [1/Pa]
    %Isothermal Compressibility of water, in units of inverse Pascals
    Beta_T_invPa_Lepple = (  5.062271                      ...
                           - 0.031790 * temperature_C      ...
                           + 0.000407 * temperature_C^2 )  ...
                           * 10^-10;
                      
    %In units of inverse bars
    %Lepple's Std Dev is 0.09 × 10?6 bar?1
    Beta_T_invBar_Lepple     = Beta_T_invPa_Lepple * 100000;
    Beta_T_invBar_Lepple_std = (0.09e-6).*ones(size(Beta_T_invBar_Lepple));

    %Convert uncertainty to [1/Pa]
    Beta_T_invPa_Lepple_std = Beta_T_invBar_Lepple_std / 100000;
    
    %In units of inverse atmospheres
    Beta_T_invAtm_Lepple     = Beta_T_invBar_Lepple * 1.01325;
    Beta_T_invAtm_Lepple_std = Beta_T_invBar_Lepple_std * 1.01325;
    
    %Use Lepples value of BetaT for calculations
    Beta_T_invPa        = Beta_T_invPa_Lepple;
    Beta_T_invPa_std    = Beta_T_invPa_Lepple_std;
    
%      ---------------------[ Method #2 ]-----------------------    
elseif(strcmpi(BetaT_Method, 'Kell'))
    %From [Kell, 1970]
    %Isothermal Compressibility in 10^6 * [1/bar] of Water
    Beta_T_Kell_invBar =  ...
                       (50.88630                                ...
                    +   0.7171582 * 10^0  * temperature_C       ...
                    +   0.7819867 * 10^-3 * temperature_C^2     ...
                    +   31.622140 * 10^-6 * temperature_C^3     ...
                    -   0.1323594 * 10^-6 * temperature_C^4     ...
                    +   0.6345750 * 10^-9 * temperature_C^5)    ...
                    /  ( 1                                      ...
                    +   21.659280 * 10^-3 * temperature_C);
    
    %Convert to [1/bar]
    %Kell Std Dev is 0.002 x 10^-6 [1/bar]
    Beta_T_Kell_invBar      = Beta_T_Kell_invBar * 10^-6;
    Beta_T_Kell_invBar_std  = (0.002e-6).*ones(size(Beta_T_Kell_invBar));
    
    %In Units of [1/Pa]
    Beta_T_Kell_invPa       = Beta_T_Kell_invBar        / 100000;
    Beta_T_Kell_invPa_std   = Beta_T_Kell_invBar_std    / 100000;
    
    %In units of [1/atm]
    Beta_T_Kell_invAtm      = Beta_T_Kell_invBar     * 1.01325;
    Beta_T_Kell_invAtm_std  = Beta_T_Kell_invBar_std * 1.01325;
    
    %Use Kell's value of BetaT for calculations
    Beta_T_invPa        = Beta_T_Kell_invPa;
    Beta_T_invPa_std    = Beta_T_Kell_invPa_std;
        
end

%==========================================================================
%%%%%%%%%% C O M P U T E  D E N S I T Y  O F  W A T E R %%%%%%%%%%%%%%%%%%%
%==========================================================================

%      ---------------------[ Method #1 ]-----------------------
if(strcmpi(Density_Method, 'BryanCox'))
%The equation for conventional density [Bryan & Cox 1972]
%Also mentioned in Mamayev, 1975, pg 45 Eqn 10.2
%sigma_T(T,S,Z) - sigma_0(T0,S0,Z) = X1*(T-T0) + X2*(S-S0) + X3*(T-T0)^2
%sigma_T, sigma_0 are conventional density = (rho - 1)10^3, rho in g/cm^3

%Depth in meters (Z)
Z_m = (0:250:6000)';

%X Coefficients as a function of depth
%               X1(Z)       X2(Z)       X3(Z)
X_coeffs = [    -0.19494    0.77475     -0.49038E-02;
                -0.15781    0.78318     -0.52669E-02;
                -0.13728    0.78650     -0.55258E-02;
                -0.12720    0.78807     -0.56610E-02;
                -0.12795    0.78710     -0.56274E-02;
                -0.12312    0.78763     -0.56972E-02;
                -0.11837    0.78822     -0.57761E-02;
                -0.11896    0.78751     -0.57631E-02;
                -0.12543    0.78560     -0.56422E-02;
                -0.13138    0.78368     -0.55239E-02;
                -0.13250    0.78300     -0.55116E-02;
                -0.13871    0.78109     -0.53946E-02;
                -0.14483    0.77920     -0.52793E-02;
                -0.15088    0.77733     -0.51654E-02;
                -0.15673    0.77544     -0.50557E-02;
                -0.15771    0.77475     -0.50466E-02;
                -0.16363    0.77292     -0.49360E-02;
                -0.16948    0.77110     -0.48268E-02;
                -0.17524    0.76930     -0.47193E-02;
                -0.18556    0.76641     -0.45102E-02;
                -0.19107    0.76467     -0.44074E-02;
                -0.19650    0.76295     -0.43061E-02;
                -0.20186    0.76126     -0.42068E-02;
                -0.20715    0.75958     -0.41138E-02;
                -0.21237    0.75792     -0.40172E-02;];

%Conventional density, mean temperature and salinity as a function of depth
%            (rho0-1)10^-3  T0          S0
sigma0T0S0  = [ 24.458      13.5        32.600;
                28.475      8.50        35.150;
                29.797      6.00        34.900;
                31.144      4.50        34.900;
                32.236      4.00        34.750;
                33.505      3.00        34.750;
                34.808      2.00        34.800;
                35.969      1.50        34.750;
                37.143      1.50        34.800;
                38.272      1.50        34.800;
                39.462      1.00        34.800;
                40.582      1.00        34.800;
                41.695      1.00        34.800;
                42.801      1.00        34.800;
                43.863      1.00        34.750;
                45.038      0.50        34.750;
                46.130      0.50        34.750;
                47.216      0.50        34.750;
                48.296      0.50        34.750;
                49.278      1.00        34.750;
                50.344      1.00        34.750;
                51.404      1.00        34.750;
                52.459      1.00        34.750;
                53.508      1.00        34.750;
                54.552      1.00        34.750;];
                
%Interpolate coefficients for the given depth
X1      = interp1(Z_m, X_coeffs(:,1),   depth_m, 'linear');
X2      = interp1(Z_m, X_coeffs(:,2),   depth_m, 'linear');
X3      = interp1(Z_m, X_coeffs(:,3),   depth_m, 'linear');
sigma0  = interp1(Z_m, sigma0T0S0(:,1), depth_m, 'linear');
T0      = interp1(Z_m, sigma0T0S0(:,2), depth_m, 'linear');
S0      = interp1(Z_m, sigma0T0S0(:,3), depth_m, 'linear');

%Formula for conventional density (rho - 1)*10^3, rho in g/cm^3
sigT = @(Tc, Sppt)(  X1*(Tc - T0)         ...
                   + X2*(Sppt - S0)       ...
                   + X3*(Tc - T0)^2       ...
                 ) + sigma0;

%Compute conventional density (rho - 1)*10^3, rho in g/cm^3
sigmaT          = sigT(temperature_C, salinity_ppt); %Density at T and S      
sigmaT_4C       = sigT(4, salinity_ppt);             %Density at T=4C and S
sigmaT_pure     = sigT(temperature_C, 0);            %Density at T and S=0    
sigmaT_pure_4C  = sigT(4, 0);                        %Density at T=4, S=0

%Convert from conventional density to density [g/cm^3]    
%sigma0 = (rho - 1)10^3 ---> rho = (sigma0/10^3)+1
h2o_density_gcm3        = (sigmaT           /10^3) + 1;
h2o_density_4C_gcm3     = (sigmaT_4C        /10^3) + 1;
pureh2o_density_gcm3    = (sigmaT_pure      /10^3) + 1;
pureh2o_density_4C_gcm3 = (sigmaT_pure_4C   /10^3) + 1;

%Convert from density in [g/cm^3] to [g/ml]
%NOTE: Sea water is more dense than pure water, so density can be greater
%than 1 [g/ml]
h2o_density_gml         = h2o_density_gcm3          * 1.000027;
h2o_density_4C_gml      = h2o_density_4C_gcm3       * 1.000027;
pureh2o_density_gml     = pureh2o_density_gcm3      * 1.000027;
pureh2o_density_4C_gml  = pureh2o_density_4C_gcm3   * 1.000027; 

%Convert to [kg/m^3]
h2o_density_kgm3        = h2o_density_gcm3        * 1000;
h20_density_4C_kgm3     = h2o_density_4C_gcm3     * 1000;
pureh2o_density_kgm3    = pureh2o_density_gcm3    * 1000;
pureh2o_density_4C_kgm3 = pureh2o_density_4C_gcm3 * 1000;
            
%      ---------------------[ Method #2 ]-----------------------
elseif(strcmpi(Density_Method, 'Mamayev'))
    
%Coefficients for Mamayev's equation (9.6)
X1 = 28.152;    X2 = 0.0735;    X3 = 0.00469;
X4 = 0.802;     X5 = 0.002;     X6 = 35;
    
%Formula for conventional density For 1 atm using Eqn 9.6
sigT = @(Tc, Sppt)     X1            ...
                    -  X2 * Tc       ...
                    -  X3 * Tc^2     ...
                    + (X4 - X5*Tc)   ...
                    * (Sppt - X6);

%Compute conventional density (rho - 1)*10^3, rho in g/cm^3
sigmaT          = sigT(temperature_C, salinity_ppt); %Density at T and S      
sigmaT_4C       = sigT(4, salinity_ppt);             %Density at T=4C and S
sigmaT_pure     = sigT(temperature_C, 0);            %Density at T and S=0    
sigmaT_pure_4C  = sigT(4, 0);                        %Density at T=4, S=0                                
  
%Convert from conventional density to density [g/cm^3]    
%sigma0 = (rho - 1)10^3 ---> rho = (sigma0/10^3)+1
h2o_density_gcm3        = (sigmaT           /10^3) + 1;
h2o_density_4C_gcm3     = (sigmaT_4C        /10^3) + 1;
pureh2o_density_gcm3    = (sigmaT_pure      /10^3) + 1;
pureh2o_density_4C_gcm3 = (sigmaT_pure_4C   /10^3) + 1;

%Convert from density in [g/cm^3] to [g/ml]
%NOTE: Sea water is more dense than pure water, so density can be greater
%than 1 [g/ml]
h2o_density_gml         = h2o_density_gcm3          * 1.000027;
h2o_density_4C_gml      = h2o_density_4C_gcm3       * 1.000027;
pureh2o_density_gml     = pureh2o_density_gcm3      * 1.000027;
pureh2o_density_4C_gml  = pureh2o_density_4C_gcm3   * 1.000027;

%Convert to [kg/m^3]
h2o_density_kgm3        = h2o_density_gcm3        * 1000;
h20_density_4C_kgm3     = h2o_density_4C_gcm3     * 1000;
pureh2o_density_kgm3    = pureh2o_density_gcm3    * 1000;
pureh2o_density_4C_kgm3 = pureh2o_density_4C_gcm3 * 1000;

%      ---------------------[ Method #3 ]-----------------------
elseif(strcmpi(Density_Method, 'UNESCO'))
  %This is the UNESCO standard method
  %This method has a Std Dev of 3.6e-3 [kg/m^3] for 1 atm
  
    %Convert depth in meters to pressure in bars
    
    %Calculate depth from a range of pressures, then interpolate to the
    %exact depth needed.  I can't find an easy equation to go from depth to
    %pressure, so I calculated depths at many pressures and interpolated to
    %the closest value.
    %This algorithm is from UNESCO
    P           = 0:0.01:10000;     %decibars
    latitude    = 42;               %No need to make this a variable
    X           = sin(latitude / 57.29578)^2;
    GR          = 9.780318 * (1.0 + (5.2788E-3 + 2.36E-5 * X) * X) ...
                + 1.092E-6 .* P;
    DepthTerm   = (((-1.82E-15 .* P + 2.279E-10) .* P - 2.2512E-5) ...
                  .* P + 9.72659) .* P;
    DEPTH       = DepthTerm ./ GR ;
    
    p_dbar = interp1(DEPTH, P, depth_m, 'linear');
    p_bar = p_dbar/10;
    
        
    %Compute standard mean ocean water density in [kg/m^3]
    %Coefficients
    w0 = 999.842594;    w1 = 6.793952e-2;    w2 = 9.095290e-3;
    w3 = 1.001685e-4;   w4 = 1.120083e-6;    w5 = 6.536332e-9;
    
    %Function for standard mean ocean water density in [kg/m^3]
    %UNESCO, 1981
    SMOW      = @(Tc)   w0                      ...
                    +   w1 * Tc      ...
                    -   w2 * Tc^2    ...
                    +   w3 * Tc^3    ...
                    -   w4 * Tc^4    ...
                    +   w5 * Tc^5;
    
    %Compute the density of seawater at 1 atm
    %Coefficients
    p0 =  8.24493e-1;    p1 =  4.0899e-3;    p2 =  7.6438e-5;
    p3 =  8.2467e-7;     p4 =  5.3875e-9;    p5 = -5.72466e-3;
    p6 =  1.0227e-4;     p7 =  1.6546e-6;    p9 =  4.8314e-4;
    
    %Formula for density of ocean water in [kg/m^3]
    h2o_density = @(Tc,Sppt)     SMOW(Tc)          ...
                            +  ( p0                ...
                               - p1 * Tc           ...
                               + p2 * Tc^2         ...
                               + p3 * Tc^3         ...
                               + p4 * Tc^4         ...
                               )    * Sppt         ...
                            +  ( p5                ...
                               + p6 * Tc           ...
                               - p7 * Tc^2         ...
                               )    * Sppt^(3/2)   ...
                            +    p9 * Sppt^2;                        
    
    %Formula for correcting the density for depth (if needed)     
    %Coefficients for Kw and K_st0
    kw0 = 19652.21;    kw1 = 148.4206;        kw2 = 2.327105;
    kw3 = 1.360477e-2; kw4 = 5.155288e-5;     k0 = 54.6746;
    k1 = 0.603459;     k2 = 1.09987e-2;       k3 = 6.1670e-5;
    k4 = 7.944e-2;     k5 = 1.6483e-2;        k6 = 5.3009e-4;
    
    %Formula for Kw
    Kw = @(Tc)  kw0                ...
              + kw1 * Tc           ...
              - kw2 * Tc^2         ...
              + kw3 * Tc^3         ...
              - kw4 * Tc^4;
    
    %Formula for Kw(T,S) at p=0 (1 atm)
    K_st0 = @(Tc,Sppt)    Kw(Tc)            ...
                     +  ( k0                ...
                        - k1 * Tc           ...
                        + k2 * Tc^2         ...
                        - k3 * Tc^3         ...
                        )    * Sppt         ...
                     +  ( k4                ...
                        - k5 * Tc           ...
                        + k6 * Tc^2         ...
                        )    * Sppt^(3/2);          
          
    %Coefficients for Aw
    aw0 = 3.239908;    aw1 = 1.43713e-3;    aw2 = 1.16092e-4;  
    aw3 = 5.77905e-7;  a0 = 2.2838e-3;      a1 = 1.0981e-5;
    a2 = 1.6078e-6;    a3 = 1.91075e-4;
    
    Aw = @(Tc)  aw0                ...
              + aw1 * Tc           ...
              + aw2 * Tc^2         ...
              - aw3 * Tc^3;

    A = @(Tc, Sppt)     Aw(Tc)            ...
                   +  ( a0                ...
                      - a1 * Tc           ...
                      - a2 * Tc^2         ...              
                      )    * Sppt         ...              
                      + a3 * Sppt^(3/2);                    
          
    %Coefficients for Bw
    bw0 = 8.50935e-5;    bw1 = 6.12293e-6;    bw2 = 5.2787e-8;    
    b0 = -9.9348e-7;     b1  = 2.0816e-8;     b2  = 9.1697e-10;
    
    Bw = @(Tc)  bw0                ...
              - bw1 * Tc           ...
              + bw2 * Tc^2;    
          
    B = @(Tc, Sppt)   Bw(Tc)             ...
                 +  ( b0                 ...
                    + b1 * Tc            ...
                    + b2 * Tc^2          ...              
                    )    * Sppt;                    
          
    %Put together the secant bulk modulus
    K_stp =@(Tc, Sppt)     K_st0(Tc,Sppt)       ...
                         + A(Tc,Sppt)*p_bar     ...
                         + B(Tc,Sppt)*p_bar^2;              
    
    %Formula for correcting density for depth P in bars
    rho_depth_kgm3 = @(Tc, Sppt, P) ...
                     h2o_density(Tc, Sppt) ...
                  / (1 - (P/K_stp(Tc, Sppt)));
    
    %Calculate the density of water at the given T, S, and P
    h2o_density_kgm3 = rho_depth_kgm3(temperature_C, salinity_ppt, p_bar);
    
    %Calculate the density of water at the T=4C, S, and P
    h20_density_4C_kgm3 = rho_depth_kgm3(4, salinity_ppt, p_bar);
    
    %Calculate the density of pure ocean water at T, S=0, and P
    pureh2o_density_kgm3 = rho_depth_kgm3(temperature_C, 0, p_bar);
    
    %Calculate the density of pure ocean water at T=4, S=0, and P
    pureh2o_density_4C_kgm3 = rho_depth_kgm3(4, 0, p_bar);   
    
    %Standard Deviation
    h2o_density_kgm3_std = (3.6e-3).*ones(size(h2o_density_kgm3));
    
    %Convert to [g/cm^3]
    h2o_density_gcm3        = h2o_density_kgm3          /1000;
    h2o_density_4C_gcm3     = h20_density_4C_kgm3       /1000;
    pureh2o_density_gcm3    = pureh2o_density_kgm3      /1000;
    pureh2o_density_4C_gcm3 = pureh2o_density_4C_kgm3   /1000;
    h2o_density_gcm3_std    = h2o_density_kgm3_std      /1000;
    
    %Convert to [g/ml]
    h2o_density_gml        = h2o_density_gcm3           * 1.000027;
    h2o_density_4C_gml     = h2o_density_4C_gcm3        * 1.000027;
    pureh2o_density_gml    = pureh2o_density_gcm3       * 1.000027;
    pureh2o_density_4C_gml = pureh2o_density_4C_gcm3    * 1.000027;
    h2o_density_gml_std    = h2o_density_gcm3_std       * 1.000027;
end

%==========================================================================
%%%%%%%%% C O M P U T E  R E F R.  I N D E X  O F  W A T E R %%%%%%%%%%%%%%
%==========================================================================

%      ---------------------[ Method #1 ]-----------------------
if(strcmpi(RefrIdx_Method, 'Eisenberg'))
%From Jonasz-Fournier 2007 and Eisenberg 1965
% NOTE B(lambda) in JF-2007 are wrong!  I recomputed the values
% New Values are given below

%               Wavelength  A               B               C 
Coeff_Eis = [   404.66,     0.2119526,      0.87642,        6.5413e-05;
                435.83,     0.2105352,      0.87898,        6.5122e-05;
                447.15,     0.2100945,      0.87974,        6.4978e-05;
                471.31,     0.2092569,      0.88110,        6.4619e-05;
                486.13,     0.2088026,      0.88180,        6.4371e-05;
                501.57,     0.2083698,      0.88249,        6.4079e-05;
                546.07,     0.2073063,      0.88410,        6.3121e-05;
                576.96,     0.2066930,      0.88503,        6.2358e-05;
                587.56,     0.2065009,      0.88532,        6.2084e-05;
                589.26,     0.2064709,      0.88538,        6.2037e-05;
                656.28,     0.2054298,      0.88694,        6.0146e-05;
                667.81,     0.2052736,      0.88714,        5.9805e-05;
                706.52,     0.2047862,      0.88796,        5.8583e-05];

%Interpolate for the given wavelength
A_eis = interp1(Coeff_Eis(:,1), Coeff_Eis(:,2), ...
                wavelength_nm, 'linear', 'extrap');
B_eis = interp1(Coeff_Eis(:,1), Coeff_Eis(:,3), ...
                wavelength_nm, 'linear', 'extrap');
C_eis = interp1(Coeff_Eis(:,1), Coeff_Eis(:,4), ...
                wavelength_nm, 'linear', 'extrap');


%Constants A,B,C as a function of wavelength (Fitted)
% % A =         1.95197E-01                         ...
% %         +   8.94225E+00 * wavelength_nm.^-1     ...
% %         -   2.40736E+03 * wavelength_nm.^-2     ...
% %         +   6.20279E+05 * wavelength_nm.^-3;
% % 
% % B =         8.92632E-01                         ...
% %         +   7.36562E-01 * wavelength_nm.^-1     ...
% %         -   2.92809E+03 * wavelength_nm.^-2     ...
% %         +   0.00000E+00 * wavelength_nm.^-3;    
% %     
% % C =         5.24970E-05                         ...
% %         +   8.27938E-08 * wavelength_nm.^1     ...
% %         -   1.53429E-10 * wavelength_nm.^2     ...
% %         +   6.85867E-14 * wavelength_nm.^3;

%Compute Eisenberg Eqn for sea water
%This is f(n',L,T) = (n'^2-1)/(n'^2+2)
f_nLT_sea = A_eis .* (h2o_density_gml...
                    /pureh2o_density_4C_gml...
                 ).^B_eis .* exp(-1 .* C_eis .* temperature_C);

%Compute Eisenberg Eqn for pure water
f_nLT_pure = A_eis .* ( pureh2o_density_gml...
                    /pureh2o_density_4C_gml...
            ).^B_eis .* exp(-1 .* C_eis .* temperature_C);

%Compute the refractive indices
n_PureWater = ((1 + 2*f_nLT_pure) ./ (1 - f_nLT_pure)).^0.5;
n_SeaWater =  ((1 + 2*f_nLT_sea)  ./ (1 - f_nLT_sea )).^0.5;

%Eisenberg's data was based on a dataset at one wavelength, so it is
%difficult to gauge a Std Dev.  But on average the error is 3e-7.
[n_PureWater_std, n_SeaWater_std] = deal((3e-7).*ones(size(n_PureWater)));

%      ---------------------[ Method #2 ]-----------------------
elseif(strcmpi(RefrIdx_Method, 'QuanFry'))
%Correct Quan and Fry number for vacuum wavelength [in water]

%Compute the refractive index of air using [Ciddor, 1996]
%Coefficients 
k0 = 238.0185;          %[um^-2]
k1 = 5792105;           %[um^-2]
k2 = 57.362;            %[um^-2]
k3 = 167917;            %[um^-2]
nu = 1./wavelength_um;  %[um^-1]

%Get refr Idx of air from Ciddor, 1996
n_air = (((k1./(k0 - nu.^2)) + (k3./(k2 - nu.^2))) / 10^8) + 1;

%Ciddor claims un uncertainty in refr idx of 2-5 x 10^-8
n_air_std = (5e-8).*ones(size(n_air));

%Compute the refractive index of water w.r.t. air using [Quan & Fry 1995]
%Coefficients
n0 = 1.31405;       n1 = 1.779e-4;      n2 = -1.05e-6;
n3 = 1.6e-8;        n4 = -2.02e-6;      n5 = 15.868;
n6 = 0.01155;       n7 = -0.00423;      n8 = -4382;
n9 = 1.1455e6;

%Quan and Fry's Refractive index formula
n_QF = @(Tc, Sppt, Wave_nm)                 ...
    n0 +                                    ...
    Sppt.*(n1 + n2.*Tc + n3.*(Tc.^2) ) +    ...
    n4.*(Tc.^2) +                           ...
    (n5 + n6.*Sppt + n7.*Tc)./Wave_nm +     ...
    n8./(Wave_nm.^2) +                      ...
    n9./(Wave_nm.^3);

%Calculate refractive index for sea water
%The data used by Quan & Fry has an uncertainty in n of 3 x 10^-5;
n_sea       = n_QF(temperature_C, salinity_ppt, wavelength_nm);
n_sea_std   = (3e-5).*ones(size(n_sea));

%Calculate refractive index for pure water
n_pure = n_QF(temperature_C, 0, wavelength_nm);
n_pure_std   = (3e-5).*ones(size(n_pure));

%Correct Quan & Frys Refr Index so it's w.r.t. vacuum.
n_PureWater  = n_pure .* n_air;
n_SeaWater   = n_sea  .* n_air;

%Compute standard deviation for n_PureWater
n_PureWater_std = sqrt( (n_pure_std./n_pure).^2     ...
                      + (n_air_std./n_air).^2)      ...
                     .* n_PureWater;
                      
%Compute standard deviation for n_SeaWater                      
n_SeaWater_std = sqrt( (n_sea_std./n_sea).^2        ...
                      + (n_air_std./n_air).^2)      ...
                     .* n_SeaWater;                      

end

%==========================================================================
%%%%%% C H A N G E  I N  R E F R.  I N D E X  W R T  P R E S S U R E %%%%%%
%==========================================================================

%Original data for dndP_663nm from O'Connor and Shulpf 1967
%           Tc      x 10^-10 [1/Pa]
Oc_Sh_67 = [05,     1.60;
            06,     1.59;
            07,     1.58;
            08,     1.57;
            09,     1.56;
            10,     1.56;
            15,     1.52;
            20,     1.50;
            25,     1.47;
            30,     1.45;
            35,     1.42;];
        
%Estimate Std Uncertainty from difference between fit and value
fitted_dndP_663nm = (1.61857 - 0.005785  .* Oc_Sh_67(:,1)) .* 10^-10;
dndP_663nm_std    = sqrt( mean( ...
                        (fitted_dndP_663nm - (Oc_Sh_67(:,2)*10^-10)).^2 ));
                                       

%From Jonasz & Fournier, 2007    %units = (1/Pa)
dndP_663nm      = (1.61857 - 0.005785  * temperature_C) * 10^-10; %Eqn 2.54
dndP_20C        = (1.5989  - 0.000156 .* wavelength_nm) * 10^-10; %Eqn 2.53
dndP_663nm_20C  = (1.5989  - 0.000156 .* 663)           * 10^-10; %Eqn 2.53
dndP_t_invPa    = (dndP_20C .* dndP_663nm) ./ dndP_663nm_20C;     %(1/Pa)

%I can't find the original text for Evtyushenkov & Kiyachenko 1982
%So I don't have an estimate of StdDev for dndP_20C. Assuming +-5% for now.
reps = 10000;
dim = firstNonSingletonDimension(dndP_20C);
if(dim==1)    
    base = dndP_20C*ones(1,reps);
    uncertainty = 0.1 .*rand(size(base)) + 0.95;% 1.05 > rand > 0.95    
    dndP_20C_std = std(base.*uncertainty, 0, 2);
else   
    base = ones(reps,1)*dndP_20C;
    uncertainty = 0.1 .*rand(size(base)) + 0.95;% 1.05 > rand > 0.95    
    dndP_20C_std = std(base.*uncertainty, 0, 1);
end

%Compute the uncertainty in dndP_t_invPa
numerator_std = sqrt( (dndP_20C_std./dndP_20C).^2 ...
                    + (dndP_663nm_std./dndP_663nm).^2 ) ...
                   .* (dndP_20C .* dndP_663nm);
dndP_t_invPa_std = ...
            sqrt( (numerator_std./mean((dndP_20C .* dndP_663nm))).^2 ...
                + (dndP_20C_std./dndP_20C).^2 ) ...
               .*  dndP_t_invPa;              

%==========================================================================
%%%%%% C H A N G E  I N  R E F R.  I N D E X  W R T  D E N S I T Y %%%%%%%%
%==========================================================================
%This is from Proutiere, A., et al. (1992). "Refractive index and density 
%variations in pure liquids: a new theoretical relation." The Journal of 
%Physical Chemistry 96(8): 3485-3489.
%Used in Zhang equation for Riso
%Uncertainties are -1.6 < % < 6.1

%For sea water
rhodn2drho_sea =  (n_SeaWater.^2 - 1)                           ...
               .* (  1                                          ...
                  +  (2/3)                                      ...
                  .* (n_SeaWater.^2 + 2)                        ...
                  .* ((n_SeaWater.^2 - 1)./(3.*n_SeaWater)).^2);

%Compute uncertainty estimate for rhodn2drho_sea
reps = 10000;
dim = firstNonSingletonDimension(rhodn2drho_sea);
if(dim==1)    
    base = rhodn2drho_sea*ones(1,reps);
    uncertainty = 0.077.*rand(size(base))+0.984;% 1.061 < rnd # < 0.984    
    rhodn2drho_sea_std = std(base.*uncertainty, 0, 2);
else   
    base = ones(reps,1)*rhodn2drho_sea;
    uncertainty = 0.077.*rand(size(base))+0.984;% 1.061 < rnd # < 0.984     
    rhodn2drho_sea_std = std(base.*uncertainty, 0, 1);
end

              
%For pure water              
rhodn2drho_pure = (n_PureWater.^2 - 1)                           ...
               .* (  1                                           ...
                  +  (2/3)                                       ...
                  .* (n_PureWater.^2 + 2)                        ...
                  .* ((n_PureWater.^2 - 1)./(3.*n_PureWater)).^2);

%Compute uncertainty estimate for rhodn2drho_pure
% reps = 10000;
% base = ones(reps,1)*rhodn2drho_pure;
% uncertainty = 0.077.*rand(size(base))+0.984;%rand # between 1.061 and 0.984
% rhodn2drho_pure_std = std(base.*uncertainty, 0, 1);
reps = 10000;
dim = firstNonSingletonDimension(rhodn2drho_pure);
if(dim==1)    
    base = rhodn2drho_pure*ones(1,reps);
    uncertainty = 0.077.*rand(size(base))+0.984;% 1.061 < rnd # < 0.984    
    rhodn2drho_pure_std = std(base.*uncertainty, 0, 2);
else   
    base = ones(reps,1)*rhodn2drho_pure;
    uncertainty = 0.077.*rand(size(base))+0.984;% 1.061 < rnd # < 0.984     
    rhodn2drho_pure_std = std(base.*uncertainty, 0, 1);
end


%==========================================================================
%%%%%% C O M P U T E  I S O T R O P I C  R A Y L E I G H  R A T I O %%%%%%%
%==========================================================================
% Riso occurrs due to minute density fluctuations of the water
% Units of Riso = [1/m]

%      ---------------------[ Method #1 ]-----------------------
if(strcmpi(RefrIdx_Method, 'Eisenberg'))
%From Jonasz & Fournier, 2007, Eqn 2.68   
%Note this equation in the book is missing a "squared" in the denominator
%            Term                                  Units
Riso =      (pi^2 ./(2* wavelength_m.^4))   ...    m^-4
        .*   Kb                             ...    Pa*m^3/K
        .*   temperature_K                  ...    K
        .*   Beta_T_invPa                   ...    1/Pa
        .*  (3 .* B_eis .* f_nLT_sea        ...
                ./(1-f_nLT_sea).^2          ...
            ).^2;     
      
%Riso for pure water        
Riso_pure = (pi^2 ./(2* wavelength_m.^4))   ...    m^-4
        .*   Kb                             ...    Pa*m^3/K
        .*   temperature_K                  ...    K
        .*   Beta_T_invPa                   ...    1/Pa
        .*  (3 .* B_eis .* f_nLT_pure       ...
                ./(1-f_nLT_pure).^2         ...
            ).^2;   
        
%Estimate Std Dev (Not exactly sure how to propagate uncertainty for the 
%last term.
lastterm        = 3 .* B_eis .* f_nLT_sea./(1-f_nLT_sea).^2;
lastterm_std    = abs(lastterm).*sqrt(2.*(n_SeaWater_std./n_SeaWater).^2);
Riso_std        = sqrt(    (lastterm_std./lastterm).^2              ...
                        +  (Beta_T_invPa_std./Beta_T_invPa).^2)     ...
                       .*  Riso;
                   
lastterm       = 3 .* B_eis .* f_nLT_pure./(1-f_nLT_pure).^2;
lastterm_std   = abs(lastterm).*sqrt(2.*(n_PureWater_std./n_PureWater).^2);
Riso_pure_std  = sqrt(    (lastterm_std./lastterm).^2              ...
                       +  (Beta_T_invPa_std./Beta_T_invPa).^2)     ...
                      .*  Riso;                   

%      ---------------------[ Method #2 ]-----------------------
elseif(strcmpi(RefrIdx_Method, 'QuanFry'))
    
%Using the method of Morel, Jonasz & Fournier
%            Term                               Units
Riso =      (2*pi^2 ./ wavelength_m.^4)  ...    m^-4
        .*   Kb                          ...    Pa*m^3/K
        .*   temperature_K               ...    K
        .*   n_SeaWater.^2               ...    
        .*   (1/Beta_T_invPa)            ...    Pa
        .*   dndP_t_invPa.^2;            ...    (1/Pa)^-2

%Formulation using Zhang 2009 (Should be more accurate)   
Riso =      (pi^2 ./ (2*wavelength_m.^4))  ...    m^-4
        .*   Kb                            ...    Pa*m^3/K
        .*   temperature_K                 ...    K
        .*   rhodn2drho_sea.^2             ...    
        .*   Beta_T_invPa;                 ...    [1/Pa]

%Compute Std Deviation for Riso (using Zhang Formulation)
Riso_std = sqrt(   (rhodn2drho_sea_std./rhodn2drho_sea).^2  ...
                +  (Beta_T_invPa_std./Beta_T_invPa).^2)     ...
                .* Riso;
    
%Riso for pure water
%Using the method of Morel, Jonasz & Fournier
Riso_pure = (2*pi^2 ./ wavelength_m.^4)  ...    m^-4
        .*   Kb                          ...    Pa*m^3/K
        .*   temperature_K               ...    K
        .*   n_PureWater.^2              ...    
        .*   (1/Beta_T_invPa)            ...    Pa
        .*   dndP_t_invPa.^2;            ...    (1/Pa)^-2    
        
%Formulation using Zhang 2009 (Should be more accurate)      
Riso_pure = (pi^2 ./ (2*wavelength_m.^4))       ...    m^-4
        .*   Kb                                 ...    Pa*m^3/K
        .*   temperature_K                      ...    K
        .*   rhodn2drho_pure.^2                 ...    
        .*   Beta_T_invPa;                      ...    [1/Pa] 
        
%Compute Std Deviation for Riso_pure (using Zhang Formulation)
Riso_pure_std = sqrt(   (rhodn2drho_pure_std./rhodn2drho_pure).^2  ...
                +  (Beta_T_invPa_std./Beta_T_invPa).^2)     ...
                .* Riso_pure;    
end

%==========================================================================
%%%%%% C H A N G E  I N  R E F R.  I N D E X  W R T  S A L I N I T Y %%%%%%
%==========================================================================
%[Quan & Fry 1995] Coefficients (Copied from Above)
n0 = 1.31405;       n1 = 1.779e-4;      n2 = -1.05e-6;
n3 = 1.6e-8;        n4 = -2.02e-6;      n5 = 15.868;
n6 = 0.01155;       n7 = -0.00423;      n8 = -4382;
n9 = 1.1455e6;

%Compute change in refractive index wrt salinity Eqn 2.59
%This is the derivative of Quan & Fry's Eqn wrt Salinity.
%Units -> S is in ppt = g/kg, so dndS = [kg/g]               
dndS_PT =  n3*temperature_C^2 + n2*temperature_C + n1 + n6./wavelength_nm;

%The error in Quan & Fry's formula is constant and independent of 
%wavelength, Tc and S.  Therefore the error of the derivative of n is 0. 
dndS_PT_std = zeros(size(dndS_PT));

%This is the derivative of n^2 with respect to salinity, using Quan and Fry
%formula. Used below in Zhang method of computing Beta_cf
dn2dS_PT = 2 * ( n3*temperature_C^2 + n2*temperature_C + n1             ...
               + n6./wavelength_nm)                                     ...
            .* ( n0                                                     ...
               + salinity_ppt                                           ...
               *   ( n3*temperature_C^2 + n2*temperature_C + n1)        ...
             + (n5 + salinity_ppt*n6 + temperature_C*n7)./wavelength_nm ...
             +  temperature_C^2*n4                                      ...
             +  n8./wavelength_nm.^2                                    ...
             +  n9./wavelength_nm.^3); 
                 
%The error in the derivative of n^2 is d/dS( (n+E)^2 ) = dn2dS +- 2EdndS
%Quan & Fry error is 3e-5.
n_QF_std        = (3e-5).*ones(size(n_SeaWater));
dn2dS_PT_std    = 2.*n_QF_std.*dndS_PT;
%==========================================================================
%%% C H A N G E  I N  S O L V.  A C T I V I T Y  W R T  S A L I N I T Y %%%
%==========================================================================        
%From Zhang (2009). "Scattering by pure seawater: Effect of salinity.

%Coefficients
a0  = -5.58651e-04;     a1  =  2.40452e-07;     a2  = -3.12165e-09;
a3  =  2.40808e-11;     a4  =  1.79613e-05;     a5  = -9.94220e-08;
a6  =  2.08919e-09;     a7  = -1.39872e-11;     a8  = -2.31065e-06;
a9  = -1.37674e-09;     a10 = -1.93316e-11;

%I think the units here are [kg/g] because S is ppt ~ [g/kg]
dlnAlphadS =    1 * (  a0                        ...
                    +  a1  * temperature_C       ...
                    +  a2  * temperature_C^2     ...
                    +  a3  * temperature_C^3     ...
                    )                            ...
            + 1.5 * (  a4                        ...
                    +  a5  * temperature_C       ...
                    +  a6  * temperature_C^2     ...
                    +  a7  * temperature_C^3     ...    
                    )      * salinity_ppt ^0.5   ...
            +   2 * (  a8                        ...
                    +  a9  * temperature_C       ...
                    +  a10 * temperature_C^2     ...
                    )      * salinity_ppt;           

%Zhang 2009 states that the uncertainty in dlnAlphadS is 0.04 percent
dlnAlphadS_std = (0.04/100).*dlnAlphadS;
                
%==========================================================================
%%%%%%%%% C O M P U T E  E L E C T R O L Y T E  E F F E C T  %%%%%%%%%%%%%%
%==========================================================================

%Salt Composition in Seawater
%From Jonasz & Fournier, 2007, Table 2.4 pg 56 -> Millero 2001 pg 252
%Mol Weight of NaCl is weird.. should be 58.44, not 55.4 as in the table.
%           Fraction    #ions   M.Weight [g/mol]    Species 
seasalts = [0.7776      2       058.44;         ... Sodium Chloride
            0.1089      3       095.21;         ... Magnesium Chloride
            0.0473      2       120.36;         ... Magnesium Sulfate
            0.0360      2       136.14;         ... Calcium Sulfate
            0.0246      3       174.26;         ... Potassium Sulfate
            0.0034      2       100.09;         ... Calcium Carbonate
            0.0022      3       184.11];        ... Magnesium Bromide

%Compute the average molecular weight per ion for the mixture.  
%This saves us having to compute each salt independently and add together
AvgWeightPerIon = sum(seasalts(:,1).*seasalts(:,3)./seasalts(:,2));

%Zhang 2009 calculated this number to be 31.33, I trust his result more
%than Jonasz & Fournier
AvgWeightPerIon = 31.33;

%Calculate the additional scattering due to concentration fluctuation [1/m]

%With Assumption of ideal solution (no dlnAdS). Better to use non-ideal
%From Jonasz & Fournier, 2007, Eqn 2.56   
%Density added because it's missing in JF & Morel74.  Added per Zhang 2009
%Units are correct, [1/m]
MSfluct_cf_ideal =                          ...
            (2*pi^2 ./( wavelength_m.^4))   ...    m^-4
        .*   AvgWeightPerIon                ...    g/mol
        .*   salinity_ppt                   ...    g/kg
        .*   n_SeaWater.^2                  ...    
        .*   dndS_PT.^2                     ...    [kg/g]^2
        ./   Na                             ...    mol
        ./   h2o_density_kgm3;              ...    m^3/kg

%The better way, doesn't assume ideal solutions, takes into account density
%Units are correct, [1/m].  This is the formulation of Zhang 2009
%            Term                                  Units
MSfluct_cf =  ...
            (pi^2 ./(2 * wavelength_m.^4))  ...    m^-4
        .*   MolWeight_h20                  ...    g/mol
        .*   salinity_ppt                   ...    g/kg
        .*   dn2dS_PT.^2                    ...    [kg/g]^2
        .*   10^-6                          ...    ml/m^3
        ./   Na                             ...    mol
        ./   h2o_density_gml                ...    ml/g 
        ./   -dlnAlphadS;                   ...    kg/g
%Note: the 10^-6 is not needed if the density is expressed in [g/m^3]

%Calculate the standard deviation of the electrolyte fluctuations
numerator =             (pi^2 ./(2 * wavelength_m.^4))  ...    m^-4
                    .*   MolWeight_h20                  ...    g/mol
                    .*   salinity_ppt                   ...    g/kg
                    .*   dn2dS_PT.^2                    ...    [kg/g]^2
                    .*   10^-6;                         ...    ml/m^3        
numerator_std = sqrt( 2.*(dn2dS_PT_std./dn2dS_PT).^2 ).*abs(numerator);
denominator =            Na                             ...    mol
                    .*   h2o_density_gml                ...    ml/g 
                    .*  -dlnAlphadS;                    ...    kg/g
denominator_std = sqrt( (h2o_density_gml_std./h2o_density_gml).^2   ...
                    +   (dlnAlphadS_std./dlnAlphadS).^2)            ...
                    .*   abs(denominator);

MSfluct_cf_std = sqrt( (numerator_std./numerator).^2                ...
                    +  (denominator_std./denominator).^2)           ...
                   .*   MSfluct_cf;
               
%If there is no salinity, the error due to salinity is 0               
MSfluct_cf_std(isnan(MSfluct_cf_std)) = 0;               


%==========================================================================
%%%%%% C O M P U T E  S C A T T E R I N G  C O E F F I C I E N T  %%%%%%%%%
%==========================================================================
        
%Compute Rayleigh Ratio due to density fluctuations of seawater [1/m]
R_df        =  Riso      .*   (6 + 6*delta)         ...
                         ./   (6 - 7*delta);
R_df_std    =  Riso_std  .*   (6 + 6*delta)         ...
                         ./   (6 - 7*delta);          

%Compute Rayleigh Ratio due to density fluctuations of purewater [1/m]          
R_df_pure           =  Riso_pure        .*   (6 + 6*delta)         ...
                                        ./   (6 - 7*delta);
R_df_pure_std       =  Riso_pure_std    .*   (6 + 6*delta)         ...
                                        ./   (6 - 7*delta);                    

%Compute additional Rayleigh Ratio due to concentration fluctuations of
%electrolyte [1/m]
R_cf        =  MSfluct_cf       .*   (6 + 6*delta)         ...
                                ./   (6 - 7*delta);
R_cf_std    =  MSfluct_cf_std   .*   (6 + 6*delta)         ...
                                ./   (6 - 7*delta);                

%Computer the total Rayleigh Ratio for just density fluctuations, and then 
%density + concentation fluctuations.  All units = [1/m]
Rpure           = R_df_pure;
Rsolution       = R_df + R_cf;
Rsolution_std   = sqrt(   (R_df_std).^2 + (R_cf_std).^2    );


%Compute the excess scattering contribution from salt concentrations.
%CHECK: Should be ~31% for typical salinity conditions
excessScattPercent      = Rsolution./Rpure - 1;
excessScattPercent_std  = sqrt( (Rsolution_std./Rsolution).^2   ...
                               + (R_df_pure_std./Rpure).^2)     ...
                               .* excessScattPercent;

%Compute scattering coefficient according to the formula in Morel, 74 Eq 11
b_pure          = (8*pi/3) .* ( R_df_pure ) .* (2+delta) ./ (1+delta);
b_pure_std      = (8*pi/3) .* R_df_pure_std .* (2+delta) ./ (1+delta);
b_seawater      = (8*pi/3) .* (R_df + R_cf) .* (2+delta) ./ (1+delta);
b_seawater_std  = (8*pi/3) .* Rsolution_std .* (2+delta) ./ (1+delta);


%==========================================================================
%%%%%%%%%%%%%% P U R E  W A T E R  A B S O R P T I O N  %%%%%%%%%%%%%%%%%%%
%==========================================================================

%Pope and Fry's data was taken at 22C
popefry_tempC = 22;

%Smith & Baker aw   Wave        aw[1/m]
SmithBaker_UV = [   200,        3.07;
                    210,        1.99;
                    220,        1.31;
                    230,        0.927;
                    240,        0.720;
                    250,        0.559;
                    260,        0.457;
                    270,        0.373;
                    280,        0.288;
                    290,        0.215;
                    300,        0.141;
                    310,        0.105;
                    320,        0.0844;
                    330,        0.0678;
                    340,        0.0561;
                    350,        0.0463;
                    360,        0.0379;
                    370,        0.0300;];
                
%Calculate Standard Deviation for Smith and Baker UV absorption data
%They claim "accuracy" to -5 < % < 25 between 300-480nm. No info to 200-300
%given. Assuming its the same.  
reps = 100000;
base = SmithBaker_UV(:,2)*ones(1,reps);
uncertainty = 0.3 .*rand(size(base)) + 0.95;%random # between 1.25 and 0.95
SmithBaker_UV_std = std(base.*uncertainty, 0, 2);
    

%Pope and Fry Absorption data
%                   Wave    aw[1/m]     StdDev          %Error
PopeFry_VIS = [     380,    0.01137,    0.0016,         14;
                    382.5,  0.01044,    0.0015,         15;
                    385,    0.009410,   0.0011,         13;
                    387.5,  0.009170,   0.0014,         16;
                    390,    0.008510,   0.0012,         15;
                    392.5,  0.008290,   0.0011,         14;
                    395,    0.008130,   0.0010,         13;
                    397.5,  0.007750,   0.0011,         15;
                    400,    0.006630,   0.00070,        11;
                    402.5,  0.005790,   0.00070,        12;
                    405,    0.0053,     0.00070,        14;
                    407.5,  0.005030,   0.00060,        13;
                    410,    0.004730,   0.00060,        13;
                    412.5,  0.004520,   0.00050,        13;
                    415,    0.004440,   0.00060,        13;
                    417.5,  0.004420,   0.00060,        14;
                    420,    0.004540,   0.00060,        14;
                    422.5,  0.004740,   0.00060,        13;
                    425,    0.004780,   0.00060,        14;
                    427.5,  0.004820,   0.00060,        13;
                    430,    0.004950,   0.00060,        12;
                    432.5,  0.005040,   0.00050,        11;
                    435,    0.0053,     0.00050,        11;
                    437.5,  0.0058,     0.00050,        10;
                    440,    0.006350,   0.00050,        9;
                    442.5,  0.006960,   0.00050,        9;
                    445,    0.007510,   0.00060,        8;
                    447.5,  0.0083,     0.00050,        7;
                    450,    0.009220,   0.00050,        6;
                    452.5,  0.009690,   0.00040,        6;
                    455,    0.009620,   0.00040,        5;
                    457.5,  0.009570,   0.00040,        5;
                    460,    0.009790,   0.00050,        6;
                    462.5,  0.01005,    0.00050,        6;
                    465,    0.01011,    0.00060,        7;
                    467.5,  0.01020,    0.00060,        6;
                    470,    0.01060,    0.00050,        6;
                    472.5,  0.01090,    0.00080,        8;
                    475,    0.01140,    0.00070,        7;
                    477.5,  0.01210,    0.00080,        8;
                    480,    0.01270,    0.00080,        7;
                    482.5,  0.01310,    0.00080,        7;
                    485,    0.01360,    0.00070,        6;
                    487.5,  0.01440,    0.00070,        6;
                    490,    0.015,      0.00070,        5;
                    492.5,  0.01620,    0.0014,         9;
                    495,    0.01730,    0.0010,         6;
                    497.5,  0.01910,    0.0014,         8;
                    500,    0.02040,    0.0011,         6;
                    502.5,  0.02280,    0.0012,         6;
                    505,    0.02560,    0.0013,         6;
                    507.5,  0.028,      0.0010,         5;
                    510,    0.03250,    0.0011,         4;
                    512.5,  0.03720,    0.0012,         4;
                    515,    0.03960,    0.0012,         4;
                    517.5,  0.03990,    0.0015,         5;
                    520,    0.04090,    0.00090,        3;
                    522.5,  0.04160,    0.0014,         4;
                    525,    0.04170,    0.0010,         4;
                    527.5,  0.04280,    0.0017,         5;
                    530,    0.04340,    0.0011,         4;
                    532.5,  0.04470,    0.0017,         5;
                    535,    0.04520,    0.0012,         4;
                    537.5,  0.04660,    0.0015,         4;
                    540,    0.04740,    0.0010,         3;
                    542.5,  0.04890,    0.0016,         4;
                    545,    0.05110,    0.0011,         3;
                    547.5,  0.05370,    0.0016,         4;
                    550,    0.05650,    0.0011,         3;
                    552.5,  0.05930,    0.0012,         3;
                    555,    0.05960,    0.0012,         3;
                    557.5,  0.06060,    0.0014,         4;
                    560,    0.06190,    0.0010,         3;
                    562.5,  0.064,      0.0015,         4;
                    565,    0.06420,    0.00090,        3;
                    567.5,  0.06720,    0.0014,         3;
                    570,    0.06950,    0.0011,         3;
                    572.5,  0.07330,    0.0017,         4;
                    575,    0.07720,    0.0011,         3;
                    577.5,  0.08360,    0.0016,         3;
                    580,    0.08960,    0.0012,         3;
                    582.5,  0.09890,    0.0016,         3;
                    585,    0.11,       0.0012,         3;
                    587.5,  0.1220,     0.0018,         3;
                    590,    0.1351,     0.0012,         3;
                    592.5,  0.1516,     0.0017,         3;
                    595,    0.1672,     0.0014,         3;
                    597.5,  0.1925,     0.0019,         3;
                    600,    0.2224,     0.0017,         3;
                    602.5,  0.2470,     0.0023,         3;
                    605,    0.2577,     0.0019,         3;
                    607.5,  0.2629,     0.0028,         3;
                    610,    0.2644,     0.0019,         3;
                    612.5,  0.2665,     0.0023,         3;
                    615,    0.2678,     0.0019,         3;
                    617.5,  0.2707,     0.0026,         3;
                    620,    0.2755,     0.0025,         3;
                    622.5,  0.2810,     0.0039,         3;
                    625,    0.2834,     0.0028,         3;
                    627.5,  0.2904,     0.0039,         3;
                    630,    0.2916,     0.0027,         3;
                    632.5,  0.2995,     0.0038,         3;
                    635,    0.3012,     0.0028,         3;
                    637.5,  0.3077,     0.0049,         3;
                    640,    0.3108,     0.0028,         3;
                    642.5,  0.3220,     0.0050,         3;
                    645,    0.3250,     0.0030,         3;
                    647.5,  0.3350,     0.0040,         3;
                    650,    0.34,       0.0030,         3;
                    652.5,  0.3580,     0.0060,         3;
                    655,    0.3710,     0.0030,         3;
                    657.5,  0.3930,     0.0060,         3;
                    660,    0.41,       0.0040,         3;
                    662.5,  0.4240,     0.0050,         3;
                    665,    0.4290,     0.0040,         3;
                    667.5,  0.4360,     0.0050,         3;
                    670,    0.4390,     0.0040,         3;
                    672.5,  0.4480,     0.0070,         3;
                    675,    0.4480,     0.0040,         3;
                    677.5,  0.4610,     0.0060,         3;
                    680,    0.4650,     0.0040,         3;
                    682.5,  0.4780,     0.0060,         3;
                    685,    0.4860,     0.0040,         3;
                    687.5,  0.5020,     0.0060,         3;
                    690,    0.5160,     0.0040,         3;
                    692.5,  0.5380,     0.0070,         3;
                    695,    0.5590,     0.0050,         3;
                    697.5,  0.5920,     0.0080,         3;
                    700,    0.6240,     0.0060,         3;
                    702.5,  0.6630,     0.0080,         3;
                    705,    0.7040,     0.0060,         3;
                    707.5,  0.7560,     0.0090,         3;
                    710,    0.8270,     0.0070,         3;
                    712.5,  0.9140,     0.011,          3;
                    715,    1.007,      0.0090,         3;
                    717.5,  1.119,      0.014,          3;
                    720,    1.231,      0.011,          3;
                    722.5,  1.356,      0.0080,         3;
                    725,    1.489,      0.0060,         3;
                    727.5,  1.678,      0.0070,         3];

%Smith and Baker NIR accurancy -15 < % < 10                
%                   Wave        aw[1/m]                
SmithBaker_NIR  = [ 730,        1.799;
                    740,        2.38;
                    750,        2.47;
                    760,        2.55;
                    770,        2.51;
                    780,        2.36;
                    790,        2.16;
                    800,        2.07;];
                
%Calculate Standard Deviation for Smith and Baker NIR absorption data
%They claim "accuracy" to -15 < % < 10 between 480-800nm. 
reps = 100000;
base = SmithBaker_NIR(:,2)*ones(1,reps);
uncertainty = 0.25 .*rand(size(base)) + 0.85;%random # between 1.1 and 0.85
SmithBaker_NIR_std = std(base.*uncertainty, 0, 2);                
                
%Hale and Querry (1973) Absorption Data
%                 Wave     aw[1/m]
HaleQuerry_IR = [825,	2.772217517;
                  850,	4.331701871;
                  875,	5.615372469;
                  900,	6.785840132;
                  925,	14.40038146;
                  950,	38.75733253];

%Standard deviation for Hale and Querry
reps = 100000;
base = HaleQuerry_IR(:,2)*ones(1,reps);
uncertainty = 0.25 .*rand(size(base)) + 0.85;%random # between 1.1 and 0.85
HaleQuerry_IR_std = std(base.*uncertainty, 0, 2);

% %Put aw data together
% aw_waves    = cat(1, SmithBaker_UV(:,1),    ...
%                      PopeFry_VIS(:,1),      ...
%                      SmithBaker_NIR(:,1),   ...
%                      HaleQuerry_IR(:,1));
% aw_dataset  = cat(1, SmithBaker_UV(:,2),    ...
%                      PopeFry_VIS(:,2),      ...
%                      SmithBaker_NIR(:,2),   ...
%                      HaleQuerry_IR(:,2));
% aw_std_dataset  = cat(1, SmithBaker_UV_std,    ...
%                          PopeFry_VIS(:,3),      ...
%                          SmithBaker_NIR_std,    ...
%                          HaleQuerry_IR_std);                 
                

% Pure Water IOPs
PureIOPs = [200,	3.07;
            205,	2.4807;
            210,	1.99;
            215,	1.6064;
            220,	1.31;
            225,	1.0894;
            230,	0.927;
            235,	0.8118;
            240,	0.72;
            245,	0.63421;
            250,	0.559;
            255,	0.5024;
            260,	0.457;
            265,	0.41501;
            270,	0.373;
            275,	0.32958;
            280,	0.288;
            285,	0.25127;
            290,	0.215;
            295,	0.17571;
            300,	0.141;
            305,	0.11895;
            310,	0.105;
            315,	0.094044;
            320,	0.0844;
            325,	0.075566;
            330,	0.0678;
            335,	0.061507;
            340,	0.0561;
            345,	0.051046;
            350,	0.0463;
            355,	0.041948;
            360,	0.0379;
            365,	0.034034;
            370,	0.03;
            375,	0.025638;
            380,	0.022;
            385,	0.020115;
            390,	0.0191;
            395,	0.018046;
            400,	0.0171;
            405,	0.016555;
            410,	0.0162;
            415,	0.015793;
            420,	0.0153;
            425,	0.01476;
            430,	0.0144;
            435,	0.014411;
            440,	0.0145;
            445,	0.014412;
            450,	0.0145;
            455,	0.015085;
            460,	0.0156;
            465,	0.015539;
            470,	0.0156;
            475,	0.016486;
            480,	0.0176;
            485,	0.018387;
            490,	0.0196;
            495,	0.022116;
            500,	0.0257;
            505,	0.03016;
            510,	0.0357;
            515,	0.042306;
            520,	0.0477;
            525,	0.049801;
            530,	0.0507;
            535,	0.052691;
            540,	0.0558;
            545,	0.059701;
            550,	0.0638;
            555,	0.067561;
            560,	0.0708;
            565,	0.073905;
            570,	0.0799;
            575,	0.091857;
            580,	0.108;
            585,	0.12772;
            590,	0.157;
            595,	0.20069;
            600,	0.244;
            605,	0.272;
            610,	0.289;
            615,	0.3007;
            620,	0.309;
            625,	0.31454;
            630,	0.319;
            635,	0.32378;
            640,	0.329;
            645,	0.3356;
            650,	0.349;
            655,	0.37371;
            660,	0.4;
            665,	0.41799;
            670,	0.43;
            675,	0.43973;
            680,	0.45;
            685,	0.46553;
            690,	0.5;
            695,	0.56641;
            700,	0.65;
            705,	0.73787;
            710,	0.839;
            715,	0.97238;
            720,	1.169;
            725,	1.4584;
            730,	1.799;
            735,	2.1354;
            740,	2.38;
            745,	2.4589;
            750,	2.47;
            755,	2.5095;
            760,	2.55;
            765,	2.5479;
            770,	2.51;
            775,	2.4463;
            780,	2.36;
            785,	2.2547;
            790,	2.16;
            795,	2.1038;
            800,	2.07;
            805,	2.0494;
            810,	2.0916;
            815,	2.2463;
            820,	2.4669;
            825,	2.819;
            830,	3.0975;
            835,	3.3423;
            840,	3.7098;
            845,	3.9835;
            850,	4.3758;
            855,	4.6308;
            860,	4.9393;
            865,	5.1535;
            870,	5.3659;
            875,	5.6114;
            880,	5.8288;
            885,	6.0074;
            890,	6.2567;
            895,	6.4708;
            900,	6.8198;
            905,	7.1028;
            910,	7.897;
            915,	9.4786;
            920,	11.171;
            925,	14.658;
            930,	19.27;
            935,	23.414;
            940,	29.309;
            945,	34.684;
            950,	38.304;
            955,	41.981;
            960,	44.081;
            965,	44.89;
            970,	45.31;
            975,	44.852;
            980,	43.71;
            985,	42.409;
            990,	41.418;
            995,	39.679;
            1000,	37.697];

%Standard deviation for Pure Water IOPs
reps = 100000;
base = PureIOPs(:,2)*ones(1,reps);
uncertainty = 0.25 .*rand(size(base)) + 0.85;%random # between 1.1 and 0.85
PureIOPs_std = std(base.*uncertainty, 0, 2);

%aw data
aw_waves = PureIOPs(:,1);

aw_dataset = PureIOPs(:,2);

aw_std_dataset = PureIOPs_std;

%Interpolate aw to the wavelengths given
a_pureh20 = interp1(aw_waves, ...
                    aw_dataset, wavelength_nm, 'linear');
                
%aw Standard Deviation
a_pureh20_std = interp1(aw_waves, ...
                        aw_std_dataset, wavelength_nm, 'linear');
%==========================================================================
%%%%%%%%% C O R R E C T  A W  F O R  T E M P  &  S A L I N I T Y  %%%%%%%%%
%==========================================================================                
%From Pegau,1997                

%Correction formula: aw_TS = aw + Psi_T*(T-Tr) + Psi_S*S
%Measured values from Pegau
%          Wave        Psi_T_pure      StdDev          Psi_T_sea   StdDev
Psi_T_meas = [...
           412,        0.0001,         0.0003,         0.0003,     0.0003;
           440,        0.0000,         0.0002,         0.0002,     0.0002;
           488,        0.0000,         0.0002,         0.0001,     0.0002;
           510,        0.0002,         0.0001,         0.0003,     0.0001;
           520,        0.0001,         0.0002,         0.0002,     0.0002;
           532,        0.0001,         0.0002,         0.0001,     0.0002;
           555,        0.0001,         0.0001,         0.0002,     0.0002;
           560,        0.0000,         0.0002,         0.0000,     0.0002;
           650,       -0.0001,         0.0001,        -0.0001,     0.0001;
           676,       -0.0001,         0.0001,        -0.0001,     0.0002;
           715,        0.0029,         0.0001,         0.0027,     0.0001;
           750,        0.0107,         0.0003,         0.0106,     0.0005;
           850,       -0.0065,         0.0001,        -0.0068,     0.0001;
           900,       -0.0088,         0.0001,        -0.0090,     0.0002;
           975,        0.2272,         0.0028,         0.2273,     0.0009];
       
%These are the Gaussian fitted values
%                       Wavelen   Psi_T
Psi_T_table = [        [(400:5:495)',zeros(20,1)];
                        500,    0.0001;
                        505,    0.0001;
                        510,    0.0002;
                        515,    0.0002;
                        520,    0.0002;
                        525,    0.0002;
                        530,    0.0001;
                        535,    0.0001;
                        540,    0.0001;
                        545,    0.0001;
                        550,    0.0001;
                        555,    0.0001;
                        560,    0.0001;
                        565,    0.0002;
                        570,    0.0002;
                        575,    0.0002;
                        580,    0.0003;
                        585,    0.0004;
                        590,    0.0006;
                        595,    0.0008;
                        600,    0.0010;
                        605,    0.0011;
                        610,    0.0011;
                        615,    0.0010;
                        620,    0.0008;
                        625,    0.0005;
                        630,    0.0002;
                        635,    0.00005;
                        640,   -0.0001;
                        645,    0.00005;
                        650,    0.00005;
                        655,    0.0001;
                        660,    0.0002;
                        665,    0.0002;
                        670,    0.0002;
                        675,    0.0001;
                        680,    0.00005;
                        685,   -0.0001;
                        690,   -0.0002;
                        695,   -0.0001;
                        700,    0.0002;
                        705,    0.0007;
                        710,    0.0016;
                        715,    0.0029;
                        720,    0.0045;
                        725,    0.0065;
                        730,    0.0087;
                        735,    0.0108;
                        740,    0.0122;
                        745,    0.0120;
                        750,    0.0106];
                    
%Values for pure water which have been fit to gaussians for in-between 
%wavelengths.       
Psi_T_wav = (400:5:750)';
Psi_T_pure_table = ...
          [zeros(size((400:5:495)'));
           1;            1;            2;            2;            2; 
           2;            1;            1;            1;            1; 
           1;            1;            1;            2;            2; 
           2;            3;            4;            6;            8; 
          10;           11;           11;           10;            8; 
           5;            2;            0;           -1;            0; 
           0;            1;            2;            2;            2; 
           1;            0;           -1;           -2;           -1; 
           2;            7;           16;           29;           45; 
          65;           87;          108;          122;          120;  
          106].*10^-4;
%Check: [TempDep_wav,    Psi_T_table]       


%Salinity Dependence    Wave     Psi_S          Std Dev
Psi_S_table = [         412,     0.00012,       0.00005;
                        440,    -0.00002,       0.00002;
                        488,    -0.00002,       0.00002;
                        510,    -0.00002,       0.00003;
                        532,    -0.00003,       0.00006;
                        555,    -0.00003,       0.00003;
                        650,     0.00000,       0.00003;
                        676,    -0.00002,       0.00002;
                        715,    -0.00027,       0.00006;
                        750,     0.00064,       0.00003;];



%Interpolate the slope to the wavelengths given
Psi_S       = interp1(Psi_S_table(:,1), Psi_S_table(:,2), ...
                      wavelength_nm, 'linear', 'extrap');
Psi_S_std   = interp1(Psi_S_table(:,1), Psi_S_table(:,3), ...
                      wavelength_nm, 'linear', 'extrap');

%Determine show much absorption to add based on the given salinity
additionalAbsDueToSalinity = Psi_S    .*salinity_ppt;
addAbsDueToSalinity_std    = Psi_S_std.*salinity_ppt;


%Determine the slope due to temperature
Psi_T       = interp1(Psi_T_table(:,1), Psi_T_table(:,2),   ...
                      wavelength_nm, 'linear', 'extrap');
Psi_T_std   = interp1(Psi_T_meas(:,1), Psi_T_meas(:,3),     ...
                      wavelength_nm, 'linear', 'extrap');            

%Determine how much absorption to add due to the given temperature
additionalAbsDueToTemperature = Psi_T  * (temperature_C - popefry_tempC);
addAbsDueToTemperature_std = Psi_T_std * (temperature_C - popefry_tempC);


%Adjust the final absorption due to temperature and salinity
a_seawater =    a_pureh20                       ...
              + additionalAbsDueToTemperature   ...
              + additionalAbsDueToSalinity;

%Compute the uncertainty in a_seawater          
a_seawater_std =   sqrt(    a_pureh20_std.^2                ...
                          + addAbsDueToTemperature_std.^2   ...
                          + addAbsDueToSalinity_std.^2); 

%If the correction reduces absorption below zero (unlikely, unless an
%unrealistic temperature is given), make it 0
a_seawater(a_seawater < 0) = 0;

%==========================================================================
%%%%%%%%%%%%%%% C A L C U L A T E  A T T E N U A T I O N  %%%%%%%%%%%%%%%%%
%==========================================================================          

%Compute attenuation of the water as the sum of absorption and scattering.
c_seawater      = a_seawater + b_seawater;

%Compute the uncertainty estimate for the attenuation.
c_seawater_std  = sqrt( a_seawater_std.^2 + b_seawater_std.^2);


end


%A subfunction that determines the first dimension of an array which is not
%1 (1 for column vectors, 2 for row vectors)
function dim = firstNonSingletonDimension(x)
    dim = min(find(size(x)~=1));
    if isempty(dim)
        dim = 1;
    end
end