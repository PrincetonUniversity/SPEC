function plot_spec(varargin)

% =========================================================================
%
% PLOT_SPEC - GUI tool to visualize SPEC output
%
%
% INPUT
% -----
% plot_spec can be used with either 0, 1 or 2 input argument.
%
% plot_spec(): Open the GUI tool without specifying a directory or a file.
%              By default, the directory is the current working directory.
%              The file (and directory) can be chosen from the GUI using a
%              file explorer.
% 
% plot_spec(filename) : Open the GUI tool and already specify SPEC hdf5
%              output file. By default, the directory is the current 
%              working directory
% 
% plot_spec(dir, filename) : Open the GUI tool and specify the directory
%              dir and the SPEC hdf5 output file.
%
%
% EXAMPLES
% --------
% plot_spec
% plot_spec('G2V02L1Fi.001.sp.h5')
% plot_spec('~/specruns/', 'G2V02L1Fi.001.sp.h5')
%
% 
% NOTES
% -----
% - If SPEC output file does not exist or is incomplete (no Poincare data),
%   errors will occur
% - Be sure that the path to SPEC routines has been added to MATLAB. (if
%   not, use "addpath ~/SPEC/Utilities/matlabtools/"
%
%
% Written by A.Baillod (2019)
%
% =========================================================================


% Test if a figure is open and give a new id for the GUI
listfig = findobj('type', 'figure');
id = num2str(length(listfig) +1);


% Read input variables
switch length(varargin)
    case 0
        dir = pwd;
        file = '-';
    case 1
        dir = pwd;
        file = varargin{1};
    case 2
        dir = varargin{1};
        file = varargin{2};
    otherwise
        error('Too many input arguments. Type "help plot_spec" for more informations')
end


% =========================================================================
% Open Figure and setup some tabs

set(0, 'units', 'pixels');
monitor = get(0, 'ScreenSize');
Mwidth = monitor(3);
Mheight = monitor(4);

MainFig = figure('Name', 'SPEC output', 'NumberTitle', 'off', ...
                  'OuterPosition', [0, 0, Mwidth, Mheight]);

PlotPosition = [.35 .15 .6 .75];
MainPlot = subplot('Position', PlotPosition, 'Tag', ['MainPlot_', id]);


% tab group

tabgp = uitabgroup(MainFig, 'Position', [.02, .02, .2, .96]);
tab1 = uitab(tabgp, 'Title', 'SPEC run infos');
tab2 = uitab(tabgp, 'Title', 'Plotting options');



% =========================================================================
% FIRST TAB

% File explorer to select input file and directory
SelectDir = uicontrol('Style', 'pushbutton', 'Parent', tab1, ...
                       'String', 'Browse to select output file', 'Units', 'normalized', ...
                       'Position', [.025, .94, .95, .04], 'Tag', ['SelectDirButton_', id], ...
                       'Callback', {@SelectDirButtonFct, id}, ...
                       'FontSize', 11);

% Directory input field
Dir = uicontrol('Style', 'edit', 'Parent', tab1, 'Units', 'normalized', ...
                'Position', [.025 .82 .95 .04], 'FontSize', 11, ...
                'String', dir, 'Tag', ['Dir_', id], ...
                'HorizontalAlignment', 'Left', ...
                'Callback', {@reset, id});
% And a text to describe it
DirText = uicontrol('Style', 'text', 'String', 'Directory:', 'Units', 'normalized', ...
                    'Position', [.025 .87 .5 .04], 'FontSize', 11, 'Parent', tab1, ...
                    'HorizontalAlignment', 'Left');

                
% File input field
RunName = uicontrol('Style', 'edit', 'Parent', tab1, 'Units', 'normalized', ...
                    'Position', [.025 .7 .95 .04], 'FontSize', 11, ...
                    'String', file, 'Tag', ['RunName_', id], ...
                    'Callback', {@reset, id}, ...
                    'HorizontalAlignment', 'Left');
% And a text to describe it
RunNameText = uicontrol('Style', 'text', 'String', 'Filename:', 'Units', 'normalized', ...
                        'Position', [.025 .75 .5 .04], 'FontSize', 11, 'Parent', tab1, ...
                        'HorizontalAlignment', 'Left');
                 
% Information bubble
RunInfos = uicontrol('Style', 'text', 'Units', 'normalized', 'Parent', tab1, ...
                     'Position', [.025 .2 .95 .45], ...
                     'String', sprintf('\n\n\n\n %s \n\n %s ','Select a run', 'for additional informations'), ...
                     'FontSize', 12, ...
                     'Tag', ['RunInfos_', id], 'Max', 10, ...
                     'BackgroundColor', 'w', ...
                     'HorizontalAlignment', 'Left', 'UserData', struct('data', 1));
                 
% Load run button
LoadButton = uicontrol('Style', 'pushbutton', 'Parent', tab1, ...
                       'String', 'Load run', 'Units', 'normalized', ...
                       'Position', [.025, .15, .95, .04], 'Tag', ['LoadButton_', id], ...
                       'Callback', {@LoadButtonFct, id}, ...
                       'FontSize', 11, 'UserData', struct('loaded', false)); 

                
               
        
                   
                   
                   
% =========================================================================
% SECOND TAB

% Plotting choice
PlotChoice = uicontrol('Style', 'popupmenu', 'Parent', tab2, 'Units', 'normalized', ...
                       'Position', [.025, .825, .95, .03], 'FontSize', 11, ...
                       'String', {'Poincare (default)', 'pressure', 'Toroidal flux', 'Poloidal flux', 'iota', ...
                                  'safety factor', 'modB', 'grid', ...
                                  'Surface current', 'Volume current'}, ...
                       'Callback', {@PlotChoiceFct, id}, 'Tag', ['PlotChoice_', id]);
    

% Plotting choice
PlotOption = uicontrol('Style', 'popupmenu', 'Parent', tab2, 'Units', 'normalized', ...
                       'Position', [.025, .725, .95, .03], 'FontSize', 11, ...
                       'String', {'No options'}, ...
                       'Tag', ['PlotOption_', id], ...
                       'Callback', {@PlotOptionFct, id}); 
                   
InfoButton = uicontrol('Style', 'pushbutton', 'Parent', tab2, ...
                       'String', 'Informations', 'Units', 'normalized', ...
                       'Position', [.025 .625 .95 .03], 'Tag', ['InfoButton_', id], ...
                       'Callback', {@InfoButtonFct, id}, 'FontSize', 11);
                   
% Autoplot radio button
AutoPlot = uicontrol('Style', 'radiobutton', 'String', 'Enable auto plot', 'Units', 'normalized', ...
                     'Position', [0.025, .9, .45, .03], 'FontSize', 12, ...
                     'Tag', ['AutoPlot_', id], 'Parent', tab2);
                   
                   
% Overlay radio button
Overlay = uicontrol('Style', 'radiobutton', 'String', 'Overlay', 'Units', 'normalized', ...
                     'Position', [0.525, .9, .45, .03], 'FontSize', 12, ...
                     'Tag', ['Overlay_', id], 'Parent', tab2);
                 
                 
% Display toroidal plane
DispTorPlane = uicontrol('Style', 'text', 'Parent', tab2, 'Units', 'normalized', ...
    'Position', [.025 .535 .95 .03], 'String', 'Toroidal plane 1/?', 'FontSize', 11, ...
    'Tag', ['DispTorPlane_', id], 'UserData', struct('TorPlane', 1));
% And some pushbuttons to move around
PlusTorPlane = uicontrol('Style', 'pushbutton', 'Parent', tab2, 'Units', 'normalized', ...
    'Position', [.475 .50 .45 .03], 'String', '+', 'Fontsize', 12, ...
    'Tag', ['PlusTorPlane_', id], 'Callback', {@PlusTorPlaneFct, id});
MinusTorPlane = uicontrol('Style', 'pushbutton', 'Parent', tab2, 'Units', 'normalized', ...
    'Position', [.025 .50 .45 .03], 'String', '-', 'Fontsize', 12, ...
    'Tag', ['MinusTorPlane_', id], 'Callback', {@MinusTorPlaneFct, id});
                    
                                       
% Display poloidal angle
DispTheta = uicontrol('Style', 'text', 'Parent', tab2, 'Units', 'normalized', ...
    'Position', [.025 .435 .95 .03], 'String', 'Theta = 180', 'FontSize', 11, ...
    'Tag', ['DispTheta_', id]);
% And some pushbuttons to move around
ThetaSlider = uicontrol('Style', 'slider', 'Parent', tab2, 'Units', 'normalized', ...
    'Position', [.025 .40 .95 .03], 'Value', 180, 'Min', 0, 'Max', 360, ...
    'Tag', ['ThetaSlider_', id], 'Callback', {@ThetaSliderFct, id});



% Plot button and its explanation text
Status = uicontrol('Style', 'text', 'String', 'Ready to plot', 'Units', 'normalized', ...
                   'Position', [.025 .01 .95 .15], 'FontSize', 11, 'Parent', tab2, ...
                   'Tag', ['Status_', id], 'BackgroundColor', 'w');
PlotButton = uicontrol('Style', 'pushbutton', 'Parent', tab2, ...
                       'String', 'Plot Poincare', 'Units', 'normalized', ...
                       'Position', [.025, .2, .95, .05], 'Tag', ['PlotButton_', id], ...
                       'Callback', {@PlotButtonFct, id}, ...
                       'FontSize', 11); 
               
                
                

                
             





end





% =========================================================================
%                           CALLBACK DEFINITION
% =========================================================================

function SelectDirButtonFct(src, event, id)
% Open file explorer to find SPEC output file and select it. Callback
% function for SelectDir button on tab 1
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure
%

    DirButton = findobj('Tag', ['Dir_', id]);
    FileButton = findobj('Tag', ['RunName_', id]);


    if (FileButton.String==string('-'))
        temp = 'output.sp.h5';
    else
        temp = char(FileButton.String);
    end

    [file, path] = uigetfile(temp);


    if (all(file==0) && all(path==0))
        warning('Please select a file')
    else
        DirButton.String = string(path);
        FileButton.String = string(file);
    end

    reset(src,event,id);

end


function reset(src,event,id)
% Resets options, select Poincare plot and set the flag loaded to false.
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure


    PlotChoice = findobj('Tag', ['PlotChoice_', id]);
    PlotChoice.Value = 1;

    PlotOption = findobj('Tag', ['PlotOption_', id]);
    PlotOption.Value = 1;
    PlotOption.String = 'No options';

    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    DispTorPlane.UserData.TorPlane = 1;
    DispTorPlane.String = 'Toroidal plane 1/?';

    ThetaSlider = findobj('Tag', ['ThetaSlider_', id]);
    ThetaSlider.Value = 180;

    DispTheta = findobj('Tag', ['DispTheta_', id]);
    DispTheta.String = 'Theta = 180';

    PlotButton = findobj('Tag', ['PlotButton_', id]);
    PlotButton.String = 'Plot Poincare';

    LoadButton = findobj('Tag', ['LoadButton_', id]);
    LoadButton.UserData.loaded = false;


end

function PlotChoiceFct(src,event,id)
% Define different possible options given which plot is required by the
% user. Callback of button PlotChoice
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    PlotButton = findobj('Tag', ['PlotButton_', id]);
    PlotChoice = findobj('Tag', ['PlotChoice_', id]);
    
    str_temp = PlotChoice.String(PlotChoice.Value);

    PlotButton.String = ['Plot ', str_temp{1}];
    
    tmp = pwd;
    Dir = findobj('Tag', ['Dir_', id]);
    RunName = findobj('Tag', ['RunName_', id]);
    cd(Dir.String);
    
    PlotOption = findobj('Tag', ['PlotOption_', id]);
    PlotOption.Value = 1;
    
    ChoiceList = PlotChoice.String;
    ToPlot = ChoiceList(PlotChoice.Value);
    ToPlot = ToPlot{1};    
    
    switch ToPlot
        case 'Poincare (default)'
            PlotOption.String = {'No option'};
        case 'Toroidal flux'
            PlotOption.String = {'Non cumulative', 'cumulative'};
        case 'Poloidal flux'
            PlotOption.String = {'Non cumulative', 'cumulative'};
        case 'iota'
            PlotOption.String = {'s-coordinate', 'R-coordinate', 'Toroidal flux', 'sqrt(toroidal flux)'};
        case 'safety factor'
            PlotOption.String = {'s-coordinate', 'R-coordinate', 'Toroidal flux', 'sqrt(toroidal flux)'};
        case 'grid'
            PlotOption.String = {'No option'};
        case 'pressure'
            PlotOption.String = {'No option'};
        case 'modB'
            % Load data
            try
                RunInfos = findobj('Tag', ['RunInfos_', id]);
                data = RunInfos.UserData.data;
            catch ME
               Status.String = ME.message;
               return
            end
            Nvol = data.input.physics.Nvol + data.input.physics.Lfreebound;
            NewString = {'All'};
            for ii = 1:Nvol
               NewString{ii+1} = ['Volume ', num2str(ii)]; 
            end
            PlotOption.String = NewString;
        case 'Surface current'
            PlotOption.String = {'No option'};
        case 'Volume current'
            PlotOption.String = {'Non cumulative', 'Cumulative'};
%         case 'B field'
%             PlotOption.String = {'All', 'radial component', 'poloidal component', 'toroidal component'};
    end
    
    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end     
    

    cd(tmp);
end

function InfoButtonFct(src,event,id)
% Define the information text given the plot the user requires. Callback of
% button InfoButton
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    tfigure = figure;
    t = uicontrol('Style', 'edit', 'Units', 'normalized', 'Position', [.05 .05 .9 .9], 'min', 0, 'max', 2, 'Parent', tfigure, ...
                 'FontSize', 14);
       
    PlotChoice = findobj('Tag', ['PlotChoice_', id]);
    
    ChoiceList = PlotChoice.String;
    ToPlot = ChoiceList(PlotChoice.Value);
    ToPlot = ToPlot{1};  
    
    switch ToPlot
        case 'Poincare (default)'
            str = {[newline, 'Poincare plot', newline, '-------------', newline, newline, ...
                   'Generate a Poincare plot at a given toroidal plane. Navigate ', ...
                   'between toroidal planes with the "+" and "-" buttons. Press', ...
                   ' on "Enable auto plot" for automatic plotting of each plane.']};
        case 'Toroidal flux'
            str = {[newline, 'Toroidal flux', newline, '-------------', newline, newline, ...
                   'Generate a bar plot of the toroidal flux in each volume. ', ...
                   'Options are cumulative or non-cumulative. The cumulative plot will sum ', ...
                   'the toroidal flux from volume 1 to i for the ith bar; non cumulative will ', ...
                   'not sum toroidal flux between volumes.', newline, newline, ...
                   'Please note that the flux is integrated from the vector potential solution, ', ...
                   'and not read from the variable tflux. Small differences may occur due to the ',...
                   'integration resolution.']};
        case 'Poloidal flux'
            str = {[newline, 'Poloidal flux', newline, '-------------', newline, newline, ...
                   'Generate a bar plot of the poloidal flux in each volume. ', ...
                   'Options are cumulative or non-cumulative. The cumulative plot will sum ', ...
                   'the poloidal flux from volume 1 to i for the ith bar; non cumulative will ', ...
                   'not sum poloidal flux between volumes.', newline, newline, ...
                   'Please note that the flux is integrated from the vector potential solution, ', ...
                   'and not read from the variable tflux. Small differences may occur due to the ',...
                   'integration resolution.']};
        case 'iota'
            str = {[newline, 'Rotational transform plot', newline, '-------------------------', newline, newline, ...
                    'Plots the rotational transform as a function of either the volume radial coordinate s, the major radius R,', ...
                    ' the toroidal flux or the square root of the toroidal flux (~minor radius). Change the x-axis with the ', ...
                    'option menu.']};
        case 'safety factor'
            str = {[newline, 'Safety factor plot', newline, '------------------', newline, newline, ...
                    'Plots the safety factor as a function of either the volume radial coordinate s, the major radius R,', ...
                    ' the toroidal flux or the square root of the toroidal flux (~minor radius). Change the x-axis with the ', ...
                    'option menu.']};
            
        case 'grid'
            str = {[newline, 'Grid plot', newline, '---------', newline, newline, ...
                    'Plots the coordinate grid. Represents how the volume radial coordinate s is chosen']};
        case 'pressure'
            str = {[newline, 'Pressure plot', newline, '-------------', newline, newline, ...
                    'Plots the pressure steps.']};
        case 'modB'
            str = {[newline, 'Modulus of magnetic field plot', newline, '----------------------------', newline, newline, ...
                    'Plots a colormap of the magnetic field strength on a poloidal cut (constant toroidal angle). Change the ', ...
                    'toroidal angle with the button "+" and "-".']};
        case 'Surface current'
            str = {[newline, 'Surface current plot', newline, '--------------------', newline, newline, ...
                    'Produce a bar plot of the toroidal current sheet in each KAM surface. Useful for comparison with', ...
                    'toroidal current constraint.']};            
        case 'Volume current'
            str = {[newline, 'Volume current plot', newline, '-------------------', newline, newline, ...
                    'Produce a bar plot of the volume current in each volume. Useful for comparison with', ...
                    'toroidal current constraint.']}; 
%         case 'B field'
%             str = {[newline, 'Magnetic field plot', newline, '-------------------', newline, newline, ...
%                     'Plots the toroidal, poloidal and radial component of the magnetic field as a function of', ...
%                     'the distance to the magnetic axis. The toroidal and poloidal angles can be changed with',...
%                     'the "+" / "-" buttons and the slider respectively.']}; 
    end
        
    set(t, 'String', str);
end

function PlotOptionFct(src,event,id)
% Automatically plot the data if Autoplot is set. Callback of
% button PlotOption
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end  
end

function LoadButtonFct(src,event,id)
% Load the data and set the flag loaded to true. Callback of
% button LoadButton
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    tmp = pwd;
    Dir = findobj('Tag', ['Dir_', id]);
    RunName = findobj('Tag', ['RunName_', id]);
    
    cd(Dir.String);
    
    RunInfos = findobj('Tag', ['RunInfos_', id]);
    if exist(RunName.String, 'file')
        
        data = read_spec(RunName.String);
        RunInfos.UserData.data = data;

        NewStr = sprintf('\n%s\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', ...
                         'Informations:', ...
                         '-----------------------', ...
                         sprintf('Igeometry: %i', data.input.physics.Igeometry), ...
                         sprintf('Istellsym: %i', data.input.physics.Istellsym), ...
                         sprintf('Lconstraint: %i', data.input.physics.Lconstraint), ...
                         sprintf('Lfreebound: %i', data.input.physics.Lfreebound), ...
                         sprintf('Nvol: %i', data.input.physics.Nvol), ...
                         sprintf('mupftol: %0.2e', data.input.physics.mupftol), ...
                         sprintf('ForceErr: %0.2e', data.output.ForceErr));
    else
        NewStr = sprintf('\n%s\n%s\n\n%s\n', ...
                         'Informations:', ...
                         '-----------------------', ...
                         'File not found!');
        RunInfos.String = NewStr;
        return
    end
        
                 
    
    RunInfos.String = NewStr;
   
    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end  
    
    LoadButton = findobj('Tag', ['LoadButton_', id]);
    LoadButton.UserData.loaded = true;
    
    
    %reset(src,event,id)
    
    cd(tmp);
end

function MinusTorPlaneFct(src,event,id)
% Select previous toroidal plane. Callback of
% button MinusTorPlane
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure
    
    Status = findobj('Tag', ['Status_', id]);
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    TorPlane = DispTorPlane.UserData.TorPlane;
    
    tmp = pwd;
    Dir = findobj('Tag', ['Dir_', id]);
    RunName = findobj('Tag', ['RunName_', id]);
    cd(Dir.String);
    
    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       return
    end
    
    % Change TorPlane
    nz          = size(data.poincare.R,2);  % # of toroidal planes
    
    if TorPlane==1
        TorPlane = nz;
    else
        TorPlane = TorPlane - 1;
    end
    
    DispTorPlane.UserData.TorPlane = TorPlane;
    DispTorPlane.String = sprintf('Toroidal plane %i/%i', TorPlane, nz);
    
    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end  

    cd(tmp);
end

function PlusTorPlaneFct(src,event,id)
% Select next toroidal plane. Callback of
% button PlusTorPlane
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure
    
    Status = findobj('Tag', ['Status_', id]);
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    TorPlane = DispTorPlane.UserData.TorPlane;
    
    tmp = pwd;
    Dir = findobj('Tag', ['Dir_', id]);
    RunName = findobj('Tag', ['RunName_', id]);
    cd(Dir.String);
    
    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       return
    end
    
    % Change TorPlane
    nz          = size(data.poincare.R,2);  % # of toroidal planes
    
    if TorPlane==nz
        TorPlane = 1;
    else
        TorPlane = TorPlane + 1;
    end
    
    DispTorPlane.UserData.TorPlane = TorPlane;
    DispTorPlane.String = sprintf('Toroidal plane %i/%i', TorPlane, nz);
    
    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end 
    
    cd(tmp);
end

function ThetaSliderFct(src,event,id)
% Select theta angle given the slider value. Callback of
% button ThetaSlider
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    DispTheta = findobj('Tag', ['DispTheta_', id]);
    ThetaSlider = findobj('Tag', ['ThetaSlider_', id]);
    
    theta = ThetaSlider.Value;
    DispTheta.String = sprintf('Theta = %0.1f', theta);
    
    AutoPlot = findobj('Tag', ['AutoPlot_', id]);
    if AutoPlot.Value
        PlotButtonFct(src,event,id)
    end 

end

function PlotButtonFct(src,event,id)
% Plot the data. Callback function of PlotButton
%
% INPUT
% -----
% src:      required, not used
% event:    required, not used
% id:       Identificator for plot_spec figure

    LoadButton = findobj('Tag', ['LoadButton_', id]);
    Status = findobj('Tag', ['Status_', id]);
    try
        if ~LoadButton.UserData.loaded
            error('ERROR: Please load data file first (Load run button)')
        end
    catch er
       Status.String = er.message;
       return
    end

    Status.String = ['Processing . . .'];
    
    
    tmp = pwd;
    Dir = findobj('Tag', ['Dir_', id]);
    cd(Dir.String);
    
    MainPlot = findobj('Tag', ['MainPlot_', id]);
    %axes(MainPlot);
    
    RunName = findobj('Tag', ['RunName_', id]);
    PlotChoice = findobj('Tag', ['PlotChoice_', id]);
    
    ChoiceList = PlotChoice.String;                                            % 'Poincare (default)', 'iota', 'safety factor'
    ToPlot = ChoiceList(PlotChoice.Value);
    ToPlot = ToPlot{1};
    
    PlotOption = findobj('Tag', ['PlotOption_', id]);
    
    % Plot
    %axes(MainPlot);
    switch ToPlot
        case 'Poincare (default)'
            out = plot_poincare(Status, id);
        case 'Toroidal flux'
            switch PlotOption.Value
                case 1
                    out = plot_torflux(Status, false, id)
                case 2
                    out = plot_torflux(Status, true, id)
            end                
        case 'Poloidal flux'
            switch PlotOption.Value
                case 1
                    out = plot_polflux(Status, false, id)
                case 2
                    out = plot_polflux(Status, true, id)
            end                  
        case 'iota'
            switch PlotOption.Value
                case 1
                    out = plot_iota(Status, 'i', 's', id);
                case 2
                    out = plot_iota(Status, 'i', 'R', id);
                case 3
                    out = plot_iota(Status, 'i', 'f', id);
                case 4
                    out = plot_iota(Status, 'i', 'r', id);
            end     
        case 'safety factor'
            switch PlotOption.Value
                case 1
                    out = plot_iota(Status, 'q', 's', id);
                case 2
                    out = plot_iota(Status, 'q', 'R', id);
                case 3
                    out = plot_iota(Status, 'q', 'f', id);
                case 4
                    out = plot_iota(Status, 'q', 'r', id);
            end     
        case 'grid'
            out = plot_grid(Status, id);
        case 'pressure'
            out = plot_pressure(Status, id);
        case 'modB'
            out = plot_modB(Status, id);
        case 'Surface current'
            out = plot_surfI(Status, id);
        case 'Volume current'
            switch PlotOption.Value
                case 1
                    out = plot_volI(Status, id, false);
                case 2
                    out = plot_volI(Status, id, true);
            end
%         case 'B field'
%             plot_field(Status, id);
    end
    
    
    cd(tmp);
    
    if out
        Status.String = 'Done. Ready to plot';
    end
end


function out = plot_poincare(Status, id)
% Plot poincare data
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%


    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
       
    if ~isfield(data, 'poincare')
        Status.String  = 'No poincare data available for this run';
        out = false;
        return
    end
        
    % Get zeta angle
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane   = DispTorPlane.UserData.TorPlane;
    Nfp         = double(data.input.physics.Nfp);
    nz          = 4*data.input.physics.Ntor*data.input.numerics.Ndiscrete;
    zeta        = double((NTorPlane-1) * (2 * pi / nz)/ Nfp);%zeta  = (nz0-1)*(2*pi/nz)/nfp;
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    plot_spec_poincare(data, NTorPlane, Nfp, 0, newfig);
    plot_spec_kam(data, zeta / (2*pi), 0);
    
    out = true;
end

function out = plot_torflux(Status, cumulative, id)
% Plot toroidal flux data. Can plot the cumulative flux (integral of B.dS
% from r=0 to r_max) or non cumulative flux (integral of B.dS
% from r_min to r_max)
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% cumulative: bool to chose if cumulative or non cumulative plot
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%

    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
        
    DispTorPlane    = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane       = DispTorPlane.UserData.TorPlane;
    Nfp             = double(data.input.physics.Nfp);
    nz              = 4*data.input.physics.Ntor*data.input.numerics.Ndiscrete;
    zeta            = double((NTorPlane-1) * (2 * pi / nz)/ Nfp);%zeta  = (nz0-1)*(2*pi/nz)/nfp;
    
    plot_spec_torflux(data, zeta, cumulative, 2);
    out = true;
end

function out = plot_polflux(Status, cumulative, id)
% Plot poloidal flux data. Can plot the cumulative flux (integral of B.dS
% from r=0 to r_max) or non cumulative flux (integral of B.dS
% from r_min to r_max)
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% cumulative: bool to chose if cumulative or non cumulative plot
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%

    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
        
    
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane = DispTorPlane.UserData.TorPlane;
    Nfp         = double(data.input.physics.Nfp);
    nz          = 4*data.input.physics.Ntor*data.input.numerics.Ndiscrete;
    zeta        = double((NTorPlane-1) * (2 * pi / nz)/ Nfp);%zeta  = (nz0-1)*(2*pi/nz)/nfp;
    
    plot_spec_polflux(data, zeta, cumulative, 2);
    out = true;
end

function out = plot_iota(Status, iorq, xaxis, id)
% Plot iota or safety factor as a function of s, R, toroidal flux or square
% root of toroidal flux.
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% iorq:     chose either the iota or safety factor profile
% xaxis:    chose s, R, toroidal flux or square root of toroidal flux for x
%           axis.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%
    
    % Load poincare data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;        
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
    
    if ~isfield(data, 'poincare')
        Status.String  = 'No poincare data available for this run. Try plotting iota at KAM surfaces using plot_spec_iota_kam.';
        out = false;
        return
    end
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    plot_spec_iota(data, iorq, xaxis, newfig);
    out = true;
end

function out = plot_grid(Status, id)
% Plot coordinate grid
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%

   % Load poincare data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
        Status.String = ME.message;
        out = false;
        return
    end
    
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane = DispTorPlane.UserData.TorPlane;
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    plot_spec_grid(data, NTorPlane, newfig);
    out = true;
end

function out = plot_pressure(Status, id)
% Plot pressure
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
        Status.String = ME.message;
        out = false;
        return
    end
    
    plot_spec_pressure(data, newfig);
    grid on;
    out = true;
end

function out = plot_modB(Status, id)
% Plot modulus of B
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%

    % Load data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
    
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane   = DispTorPlane.UserData.TorPlane;
    Nfp         = double(data.input.physics.Nfp);
    nz          = 4*data.input.physics.Ntor*data.input.numerics.Ndiscrete;
    zeta        = double((NTorPlane-1) * (2 * pi / nz) * 1 / Nfp);
    
    Mvol = data.input.physics.Nvol + data.input.physics.Lfreebound;
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    PlotOption = findobj('Tag', ['PlotOption_', id]);
    lvolstr = PlotOption.String{PlotOption.Value};
    if strcmp(lvolstr,'All')
        tarr = linspace(0, 2*pi, 100);
        sarr = linspace(-1+1E-2, 1, 100);    



        plot_spec_modB(data, 1, sarr, tarr, zeta, newfig);

        sarr = linspace(-1, 1, 100);
        for ii=2:Mvol
            plot_spec_modB(data, ii, sarr, tarr, zeta, 0);
        end

        plot_spec_kam(data, zeta / (2*pi), 0);
    else
        lvol = str2double(lvolstr(8:end));
        
        tarr = linspace(0, 2*pi, 100);
        
        if lvol==1
            sarr = linspace(-1+1E-2, 1, 100); 
        else
            sarr = linspace(-1, 1, 100); 
        end
        
        plot_spec_modB(data, lvol, sarr, tarr, zeta, newfig);
        
    end
    out = true;
end

function out = plot_surfI(Status, id)
% Plot surface current
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%
    
%Load poincare data
    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end
    
    DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
    NTorPlane   = DispTorPlane.UserData.TorPlane;
    Nfp         = double(data.input.physics.Nfp);
    nz          = 4*data.input.physics.Ntor*data.input.numerics.Ndiscrete;
    zeta        = double((NTorPlane-1) * (2 * pi / nz) * 1 / Nfp);
    
    Overlay = findobj('Tag', ['Overlay_', id]);
    if Overlay.Value
        newfig = 0;
    else
        newfig = 2;
    end
    
    plot_spec_surfcurent(data, 50, 50, zeta, newfig);
    out = true;
end

function out = plot_volI(Status, id, cumul)
% Plot volume current, either in cumulative or non cumulative format
%
% INPUT
% -----
% Status:   object containing plotting status. string is modified depending
%           on status.
% id:       Identificator for plot_spec figure
% cumul:    Chose if the cumulative volume current or the non cumulative
%           one is plotted
%
% OUTPUT
% ------
% out:      Bool set to false if failure
%

    try
        RunInfos = findobj('Tag', ['RunInfos_', id]);
        data = RunInfos.UserData.data;
    catch ME
       Status.String = ME.message;
       out = false;
       return
    end

    plot_spec_Ivolume(data, cumul, 2);
    out = true;
end

% function plot_field(Status, id)
% %Load poincare data
%     try
%         RunInfos = findobj('Tag', ['RunInfos_', id]);
%         data = RunInfos.UserData.data;
%     catch ME
%        Status.String = ME.message;
%        return
%     end
%     
%     pdata = pdata_from_data(data);
%     
%     DispTorPlane = findobj('Tag', ['DispTorPlane_', id]);
%     NTorPlane   = DispTorPlane.UserData.TorPlane;
%     Nfp         = double(pdata.Nfp);
%     nz          = size(pdata.R_lines,2);  % # of toroidal planes
%     zeta        = double((NTorPlane-1) * (2 * pi / nz) * 1 / Nfp);
%     
%     ThetaSlider = findobj('Tag', ['ThetaSlider_', id]);
%     theta = ThetaSlider.Value;
%     
%     theta = theta * pi / 180;
%     
%     Overlay = findobj('Tag', ['Overlay_', id]);
%     if Overlay.Value
%         newfig = 0;
%     else
%         newfig = 2;
%     end
%     
%     PlotOption = findobj('Tag', ['PlotOption_', id]);
%     Pltopt = PlotOption.String{PlotOption.Value};
%     
%     switch Pltopt
%         case 'All'
%             plot_spec_Bfield(data, 'all', theta, zeta, 500, newfig)
%         case 'radial component'
%             plot_spec_Bfield(data, 'psi', theta, zeta, 500, newfig)
%         case 'poloidal component'
%             plot_spec_Bfield(data, 'theta', theta, zeta, 500, newfig)
%         case 'toroidal component'
%             plot_spec_Bfield(data, 'phi', theta, zeta, 500, newfig)
%     end
%     
% end





