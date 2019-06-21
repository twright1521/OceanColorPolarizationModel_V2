function [newRXP_file_name] = Generate_RayXP_Scattering_Matrix_File_2(wavelength_air,...
                                    Re_ref_idx,Im_ref_idx,PSD,mie_file_path,file_name)

 %Delete any existing scattering matrices from unfinished processes
 if exist('mie_sc','file')   
     delete('mie_sc');    
 end
 
 % This is a hydrosol, wavelength in air must be converted
 % to wavelength in water
 refr_idx_water = getRefractiveIndex(wavelength_air);
 wavelength_water = wavelength_air ./ refr_idx_water;

 %   SizeDistIndex               PAR1        PAR2        PAR3                        
 %  ---------------------------------------------------------                                             
 %   6 (power law)               alpha       rmin        rmax                        

 parameters = {'numParticles',   1,                                ...
               'sizeDistIndex',  6,                                ...   
               'lambda',         [wavelength_water/1000, 0, 0, 0], ...
               'Rem',            [Re_ref_idx,            0, 0, 0], ...
               'Imm',            [Im_ref_idx,            0, 0, 0], ...
               'PAR1',           [PSD,                   0, 0, 0], ...
               'PAR2',           [0.1,                   0, 0, 0], ...
               'PAR3',           [50,                    0, 0, 0], ...
               'GSF',            1};
              
 % File names and paths
 outFilename  = strcat(file_name,'_wave_',num2str(wavelength_air),'.mie');
 
 RXPFormatDir = strcat(mie_file_path, filesep, 'Sol Files');
 
 % Calculate Mie Scattering Matrix
 MieOutputFilePath  = generate1ScatMatrix_mie_2(outFilename,mie_file_path,parameters{:});
 
 %Convert it to RayXP Format
 [newRXP_file_name] = MieToRXP(MieOutputFilePath, file_name, RXPFormatDir);
 
end
 
 

 
 