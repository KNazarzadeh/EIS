clear all;
close all;

load("batteryStruct_2.mat");
impedance = batteryStruct.('battery_1').impedance;
Z = load("Z_bat_tot.mat");

% Plot Nyquist diagram
% figure();
% plot(real(impedance), -imag(impedance), '-ok' , 'LineWidth' , 2 );

colors = lines(length(fieldnames(batteryStruct))); % Generate distinct colors

for i = 1:length(fieldnames(batteryStruct))
    impedance = batteryStruct.(['battery_' int2str(i)]).impedance;

    % Convert to milliohms
    % real_impedance = 1000 * real(impedance);
    % imag_impedance = 1000 * imag(impedance);
    
    real_impedance = real(impedance);
    imag_impedance = imag(impedance);

    plot(real_impedance, -imag_impedance, 'o-', ...
         'Color', colors(i, :), ...
         'MarkerFaceColor', colors(i, :), ...    
         'LineWidth', 2);
    axis equal;
    hold on;

    % Determine label text
    if i == 1
        labelText = '1 particle';
    else
        labelText = sprintf('%d particles', i);
    end

    % Add label at the end point of the curve
    % Add small offset to move text slightly away from curve
    offsetX = 0.4;
    offsetY = 0;
    lastRealLabelX = real_impedance(1) - offsetX;
    lastRealLabelY = -imag_impedance(1) + offsetY;
  
    text(lastRealLabelX, lastRealLabelY, labelText, ...
         'Color', colors(i, :), ...
         'FontSize', 14, ...
         'FontWeight', 'bold', ...
         'VerticalAlignment', 'top', ...
         'HorizontalAlignment', 'right');
end

% Place parameter text box on plot
% Get current axis limits
xLimits = xlim;
yLimits = ylim;

% Set position relative to max limits (adjust the scaling factors as needed)
xPos = xLimits(1) + 0.03 * (xLimits(2) - xLimits(1));
yPos = yLimits(1) + 0.95 * (yLimits(2) - yLimits(1));

text(xPos, yPos, 'SOC = 50%', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
text(xPos, yPos*0.96, 'T = 25^\circC', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');

text(xPos, yPos*0.03, 'f_{max} = 1 kHz', 'FontSize', 14);
text(xPos, yPos*0.85, 'f_{min} = 200 \muHz', 'FontSize', 14);

% xlabel('Re(impedance) [m\Omega]', 'Interpreter', 'tex', 'FontSize', 14);
% ylabel('-Im(impedance) [m\Omega]', 'Interpreter', 'tex', 'FontSize', 14);
xlabel('Re(impedance)', 'Interpreter', 'tex', 'FontSize', 14);
ylabel('-Im(impedance)', 'Interpreter', 'tex', 'FontSize', 14);
title('Nyquist Plot of Full-Cell Impedance - ESPM Model', 'FontSize', 18);

% Increase tick labels font size
ax = gca;
ax.FontSize = 14;

grid on;
hold off;

disp("------------------------")
disp("Completed!")

%% Plot only 5 particles
% figure();
% plot(real_impedance, -imag_impedance, 'o-', ...
%          'Color', colors(i, :), ...
%          'MarkerFaceColor', colors(i, :), ...    
%          'LineWidth', 2);
% hold on;
% labelText = sprintf('%d particles', 5);
% 
% 
% % Add label at the end point of the curve
% % Add small offset to move text slightly away from curve
% offsetX = 0.4;
% offsetY = 0;
% lastRealLabelX = real_impedance(1) - offsetX;
% lastRealLabelY = -imag_impedance(1) + offsetY;
% 
% text(lastRealLabelX, lastRealLabelY, labelText, ...
%      'Color', colors(i, :), ...
%      'FontSize', 14, ...
%      'FontWeight', 'bold', ...
%      'VerticalAlignment', 'top', ...
%      'HorizontalAlignment', 'right');
% 
% % Place parameter text box on plot
% % Get current axis limits
% xLimits = xlim;
% yLimits = ylim;
% 
% % Set position relative to max limits (adjust the scaling factors as needed)
% xPos = xLimits(1) + 0.03 * (xLimits(2) - xLimits(1));
% yPos = yLimits(1) + 0.95 * (yLimits(2) - yLimits(1));
% 
% text(xPos, yPos, 'SOC = 50%', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
% text(xPos, yPos*0.96, 'T = 25^\circC', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
% 
% text(xPos, yPos*0.03, 'f_{max} = 1 kHz', 'FontSize', 14);
% text(xPos, yPos*0.85, 'f_{min} = 400 \muHz', 'FontSize', 14);
% 
% % xlabel('Re(impedance) [m\Omega]', 'Interpreter', 'tex', 'FontSize', 14);
% % ylabel('-Im(impedance) [m\Omega]', 'Interpreter', 'tex', 'FontSize', 14);
% xlabel('Re(impedance)', 'Interpreter', 'tex', 'FontSize', 14);
% ylabel('-Im(impedance)', 'Interpreter', 'tex', 'FontSize', 14);
% title('Nyquist Plot of Full-Cell Impedance - ESPM Model', 'FontSize', 18);
% 
% % Increase tick labels font size
% ax = gca;
% ax.FontSize = 14;
% 
% grid on;
% hold off;

disp("------------------------")
disp("Completed!")