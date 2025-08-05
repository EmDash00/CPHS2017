function hcps = hcpsTools()
hcps.csv2struct = @csv2struct;
hcps.plotstruct = @plotstruct;
end


% extract the data and put into a structure
function s = csv2struct(csvFile)
    data = csvread(csvFile, 11, 0);
    data = data.';
    s.t = data(1,:); % time 
    s.treal = data(2,:); % real time
    s.r = data(3,:); % reference
    s.u = data(4,:); % inputs
    s.d = data(5,:); % disturbances
    s.y = data(6:end,:); % output
    s.e = s.r - s.y(1,:); % error
    s.data = data; % data
    s.name = csvFile;
end

function h = plotstruct(s)
set(1,'defaulttextinterpreter','latex')
timePlotOpts = {'Xlim', [5 45], 'YLim', [-1 1]};

R_COLOR = [0.9 0.7 0.1];
Y_COLOR = [0.5 0 0.6];

if length(s)>1
    STATS = true;
else
    STATS = false;
end
    
if STATS
    [u ud] = fieldStat(s, 'u');
    [y yd] = fieldStat(s, 'y');
    [e ed] = fieldStat(s, 'e');
    [r rd] = fieldStat(s, 'r');
    [d dd] = fieldStat(s, 'd');
end

t = s(1).t;

h(1) = subplot(3,2,1)
plot(t, r, 'Color', R_COLOR); hold on;
plot(t, y, 'Color', Y_COLOR);
patch([t fliplr(t)], [y-yd fliplr(y+yd)], Y_COLOR, 'FaceAlpha', 0.4);

set(gca, timePlotOpts{:});
xlabel('tracking');
ylabel('output ($y$)');

h(2) = subplot(3,2,2)
plot(t, d, 'Color', R_COLOR); hold on;
plot(t, e, 'Color', Y_COLOR);

set(gca, timePlotOpts{:});
xlabel('tracking');
ylabel('error ($r-y$)');
subplot 323

% subplot 324
% 
% subplot 325
% 
% subplot 326

fprintf('not done')
end

function [xmean xdev] = fieldStat(S, field)
    for k = 1:length(S)
        Xstack(:,:,k) = getfield(S(k), field);
    end
    xmean = mean(Xstack, 3);
    xdev = std(Xstack, [], 3);
end