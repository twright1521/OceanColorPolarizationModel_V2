%% |||||||||||||||DISCRIPTION |||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% This program extracts data for the specified inputs from the 
% RayXP output files produced from the generated input files
% from WriteRayXPinputs.m

% Output files must be generated in RayXP and put into the folder 
% 'Output Files' within their respective RayXP folder 
% of the format 'RayXP_CHL_X_WS_Y_V_Z'

clearvars

%% |||||||||||||||INFO TO CHANGE AT THE START OF EVERY SIMULATION |||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

sim_num = 8;        %Simulation Folder Number

% Angles
% Should be found in the respective angle files used to create input files
Sun_Zen = [30,45,60];       % Sun polar angle, degrees                   
Rec_Zen = [10,20,40,60];       % Receiver polar angle, degrees (integer)
Rec_Azm = 90;       % Receiver azimuth relative to sun, degrees (integer) 

% Variables to Plot

plot_vars = {'Intensity','Q','U','V','Polarization'};

[~,n] = size(plot_vars);

%% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  P R O G R A M  S T A R T %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fignum = 1;

sim_fol_name = strcat('Simulation_',num2str(sim_num));          %Simulation Folder

B = dir(sim_fol_name);

[m,~] = size(B);

for x = 3: m              %Folder from which output files will be taken
    filename{1,x-2} = string(B(x).name);     
end

for sun_zen_idx = 1:length(Sun_Zen)
    for rec_zen_idx = 1:length(Rec_Zen)
        Sim_Data_Cell = cell(m-2,7);

        k = 0;

        for w = 1: m-2

            % Load Text File

            flag = strfind(filename{w},'RayXP');

            if (isempty(flag))
                Sim_Data_Cell(w-k,:) = [];
                k = k + 1;
                continue
            else
                FID1 = fopen(strcat(sim_fol_name,filesep,filename{w},filesep,'Test_Inputs.txt'), 'rt');

                s = textscan(FID1, '%s', 'delimiter', '\n');

                fclose(FID1);

                A = string(s{1});

                for x = length(A):-1:1

                    if (A(x) == '')
                        A(x,:) = [];
                    end
                end

                A = cellstr(A);

                % Puts constants into the Sim Data cell

                cellnum = 1;

                for y = 1:9

                    if (y == 2 || y == 4)
                        continue
                    end

                    test_str = A{y};
                    idx1 = strfind(test_str, ':');

                    for x = 1:idx1
                        test_str(1) = [];
                    end

                    if (y == 1)
                        idx2 = strfind(test_str, 'mg/m^3');
                    else
                        if (y == 5)
                            idx2 = strfind(test_str, 'm/s');
                        else 
                            if (y ==6 || y == 7 || y ==8)
                                idx2 = strfind(test_str, 'degrees');
                            end
                        end
                    end

                    if ( y == 3)
                        N = str2double(test_str);
                    else
                        if ( y == 9)
                            Sim_Data_Cell{w-k,cellnum} = str2double(test_str);
                            cellnum = cellnum + 1;
                        else 
                            for x = length(test_str):-1:idx2
                                test_str(x) = [];
                            end

                            Sim_Data_Cell{w-k,cellnum} = str2double(test_str);
                            cellnum = cellnum + 1;
                        end
                    end
                end

                Sim_Data_Cell{w-k,3} = Sun_Zen(sun_zen_idx);
                Sim_Data_Cell{w-k,4} = Rec_Zen(rec_zen_idx);
                Sim_Data_Cell{w-k,5} = Rec_Azm;

                D = cell(13,N);

                % Putting TauMol, TauSol, etc into D cell

                for y = 1:N

                    for z = 1:9

                        if (z == 2)
                            continue
                        end

                        test_str = A{11*y + z};
                        idx1 = strfind(test_str, ':');

                        for x = 1:idx1
                            test_str(1) = [];
                        end

                        if (z == 1)
                            idx2 = strfind(test_str, 'nm');

                            for x = length(test_str):-1:idx2
                                test_str(x) = [];
                            end

                            D{1,y} = str2double(test_str);
                        else

                            D{z+4,y} = str2double(test_str);
                        end
                    end
                end

                % Getting I,Q,U,V,P from output file and putting into D cell

                for y = 1:N

                    outfile = strcat(sim_fol_name,filesep,filename{w},filesep,'Output Files',filesep,sprintf('Test_Script_%u.out',D{1,y}));

                    [Results] =  SingleValueRayXP(outfile,Sun_Zen(sun_zen_idx),Rec_Azm,Rec_Zen(rec_zen_idx));

                    D{2,y} = Results(3);
                    D{3,y} = Results(4);
                    D{4,y} = Results(5);
                    D{5,y} = Results(6);
                    D{6,y} = Results(7);

                end

                % Puts D cell into Sim Data cell
                Sim_Data_Cell{w-k,7} = D;

                clearvars D
            end
        end

        Sim_Data_Cell = sortrows(Sim_Data_Cell,[1 2 3 4 5 6]);
        
        save(strcat(sim_fol_name,filesep,'Sim_Data_SZ_',num2str(Sun_Zen(sun_zen_idx)),'_RZ_',num2str(Rec_Zen(rec_zen_idx)),'.mat'),'Sim_Data_Cell')

        %% ||||||||||||||| PLOT CONSTANTS |||||||||||||||||||||||||||||||||||||||||
        %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        % Variable Names

        var_names = {'[CHL]', 'Wind Speed', 'Sun Zenith', 'Receiver Zenith',...
                    'Receiver Azimuth','SSA Ocean Atm','Wavelength', 'Intensity','Q',...
                    'U','V','Polarization' 'TauMol', 'TauSol',...
                    'TauAbs', 'SSASol', 'A', 'B', 'Bb'};

        %% ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%  P L O T T I N G  S T A R T %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        [meg,~] = size(Sim_Data_Cell);

        % Make Figures Folder
        Fig_Fol = strcat(sim_fol_name,filesep,'Figures_SZ_',num2str(Sun_Zen(sun_zen_idx)),'_RZ_',num2str(Rec_Zen(rec_zen_idx)));
         if (exist(Fig_Fol, 'dir') == 0)
                mkdir(Fig_Fol)
         end

        % Create Legend and Half of Title
        y = 1;
        z = 1;
        for x = 1:6

            for w = 1:meg
                temp_array(w) = Sim_Data_Cell{w,x};
            end

            if (temp_array == temp_array(1))
                title_array(z) = temp_array(1);
                title_name{z} = var_names{x};
                z = z + 1;
            else
                legend_array(y,:) = temp_array;
                legend_name{y} = var_names{x};
                y = y + 1;
            end

            clearvars temp_array
        end


        % Title
        title_half = 'Constants -';

        for x = 1:z-1

            if (mod(x,2) && x ~= 1)
                title_half = sprintf('%s\n%s',title_half,...
                    strcat(title_name{x},':',32,num2str(title_array(x))));
            else
                title_half = strcat(title_half,32,title_name{x},':',32,num2str(title_array(x)));
            end
        end

        % Legend

        legend_cell = cell(1,meg);

        for w = 1:meg
            
            for x = 1:y-1

                if (mod(x,2) && x ~= 1)
                    legend_cell{w} = sprintf('%s\n%s',legend_cell{w},...
                        strcat(legend_name{x},':',32,num2str(legend_array(x,w))));
                else
                    legend_cell{w} = strcat(legend_cell{w},32,legend_name{x},':',32,num2str(legend_array(x,w)));
                end
            end
        end

        clearvars y z
        %% |||||||||||||||PLOTTING ||||||||||||||||||||||||||||||||||||||||||||||||

        for x = 1:n

            title_full = sprintf('%s vs. %s\n%s',plot_vars{x},var_names{7},title_half);

            idx = find(contains(var_names, plot_vars{x}))-6;

            figure(fignum)
            hold on 
            grid on

            for w = 1:meg
                temp_array = Sim_Data_Cell{w,7};

                [~,k] = size(temp_array);

                for y = 1:k
                    wav_array(y) = temp_array{1,y};
                    var_array(y) = temp_array{idx,y};
                end

                plot(wav_array,var_array,'LineWidth',0.8)
                clearvars temparray
            end

            xlabel('Wavelength - nm')

            if (contains('Intensity',plot_vars{x}))
                ylabel(strcat(plot_vars{x},32,'- W/(m^2*um*Sr)'))
            else
                if (contains('Polarization',plot_vars{x}))
                    ylabel(strcat('Degree of',32,plot_vars{x},32,'W/(m^2*um*Sr)'))
                else
                    if (contains('Q', plot_vars{x}) || ...
                            contains('U', plot_vars{x}) || ...
                            contains('V', plot_vars{x}))
                        ylabel(strcat(plot_vars{x},32,'Parameter - W/(m^2*um*Sr)'))
                    else
                        if (contains('A', plot_vars{x}) || ...
                                contains('B', plot_vars{x}) || ...
                                contains('Bb', plot_vars{x}))
                            ylabel(strcat(plot_vars{x},32,'Coefficient'))
                        else
                            ylabel(plot_vars{x})
                        end
                    end
                end
            end

            title(title_full)
            legend(legend_cell)

            fignum = fignum + 1;

            savefig(strcat(Fig_Fol,filesep,plot_vars{x}))
            
        end
        
        clearvars -except sim_num Sun_Zen Rec_Zen Rec_Azm plot_vars m n sim_fol_name filename sun_zen_idx rec_zen_idx fignum
    end
end