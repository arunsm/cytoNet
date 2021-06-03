%Save adjacency matrix and locations to .mat file
x = cellinfo(:,1,1);
y = cellinfo(:,2,1);
% adjfile = outputFileName(filename, '.mat');
% save(adjfile, 'adj', 'x', 'y');

%Show best correlated signals and uncorrelated signals
figure('visible', 'off');
set(gcf, 'color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 2]);
plot(linspace(1,timestop,framestop), squeeze(cellinfo(bestsig1,3,:)), 'r', 'LineWidth', 2);
hold on
plot(linspace(1,timestop,framestop), squeeze(cellinfo(bestsig2,3,:)), 'b', 'LineWidth', 2);
hold off
title('Best Correlated Signals');
[filepath,name,ext] = fileparts(fileName);
day=dateCode;
title({'\fontsize{20}' 'Best Correlated Signals','\fontsize{6}' '\color[rgb]{0 0.4470 0.7410}' name, day});
xlabel('Image Frame Number');
ylabel('dF/F');
%text(1,1, filename,'HorizontalAlignment','right', 'Units', 'normalized');
ax = gca;
ax.FontSize = 18;
print(sprintf('%s%s%sbestcorr', outputDirName, filesep, name),'-dtiffn')
figure('visible', 'off');
set(gcf, 'color', 'w');
plot(linspace(1,timestop,framestop), squeeze(cellinfo(badsig1,3,:)), 'r', 'LineWidth', 2);
hold on
plot(linspace(1,timestop,framestop), squeeze(cellinfo(badsig2,3,:)), 'b', 'LineWidth', 2);
% title('Bad Signals');
%xlabel('Time(min)');
xlabel('Image Frame Number');
ylabel('dF/F');
%print(outputFileName(filename,sprintf('Bad-Correlation-%i',timestop)), '-dpng')
title({'\fontsize{20}' 'Poorly Correlated Signals','\fontsize{6}' '\color[rgb]{0 0.4470 0.7410}' name, day});
ax = gca;
ax.FontSize = 18;
print(sprintf('%s%s%spoorcorr', outputDirName, filesep, name),'-dtiffn')

%Genearate heatmap from adjacency matrix
% figure, imagesc(adj);
% colorbar
% print(outputFileName(filename, sprintf('Correlation Matrix-%i',timestop)), '-dpng')

%Remove indices lower than cutoff
indices = find(adj < cutoff);
adj(indices) = 0;

%Calculate mean correlation of all cell combinations, including inactive
%ones
%inactivecombos = nchoosek(CC.NumObjects, 2) - length(corrs);
%allcorrs = [squeeze(corrs) squeeze(zeros(inactivecombos, 1))];
mncorr = mean(corrs, 'omitnan');

%Generate graph from adjacency matrix
%figure(7), imshow(segmentOverlay(bw, image), 'Border', 'Tight');
%title(sprintf('Activity = %4.3f', totalactivity/CC.NumObjects));
%figure(7),print(outputFileName(filename,'cells'), '-dpng');

[B, ~] = bwboundaries(bw, 'noholes');

figure('visible', 'off');
imshow(imadjust(image), 'Border', 'Tight');
hold on; 

for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:, 2), boundary(:, 1), 'r','LineWidth',1);
end

set(findall(gcf,'type','line'), 'LineWidth', 1.5);

[~, ~] = wgPlot(adj, [cellinfo(:,1,1) cellinfo(:,2,1)]);
title({'\fontsize{6}' '\color[rgb]{0 0.4470 0.7410}' name, day});

print(sprintf('%s%s%snetwork', outputDirName, filesep, name),'-dtiffn')

%gplot(adj, [cellinfo(:,1,1) cellinfo(:,2,1)], 'y');
%figure(7),hold on, scatter(cellinfo(:,1,1), cellinfo(:,2,1), '*g');
%hold off
%title(sprintf('Mean Correlation = %4.3f', mncorr));
%figure(7), print(outputFileName(filename,sprintf('graph-%i',timestop)), '-dpng')

%Plot correlation as a function of distance
coord = [cellinfo(:,1,1) cellinfo(:,2,1)];
%title1=sprintf('Mean Correlation = %4.3f', coord);
distInPixels = pdist(coord, 'Euclidean');
%micronsPerPixel = 0.33;
micronsPerPixel = 1;
distInMicrons = micronsPerPixel*distInPixels;

figure('visible', 'off');
set(gcf, 'color', 'w');
plot(distInMicrons, ccov, 'k.', 'MarkerSize', 20);
r = refline(0,cutoff);
r.Color = 'r';
r.LineWidth = 2;

ax = gca;
ax.FontSize = 18;
ax.YLim = [0 1];
title({'Correlation Coefficient vs. Distance','\fontsize{6}' '\color[rgb]{0 0.4470 0.7410}' name, day});
xlabel('Distance (\mum)');
ylabel('Correlation Coefficient');
print(sprintf('%s%s%scorrelation', outputDirName, filesep, name),'-dtiffn')
