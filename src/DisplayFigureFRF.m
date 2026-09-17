function DisplayFigureFRF(ExpFRF,ToolCoupledRCSAFRF,ToolRestrainedFRF,FigOpts)
% Display results

%% Variable Information


%% Experimental FRF
% Widescreen PowerPoint slide: 13 1/3 x 7.5 inches (16:9).
% Match the on-screen aspect ratio and the paper size used for JPG saving.
SlideSize = [40/3 7.5]; % ppt slide side 16:9
figure('Color','w', ...
    'Units','pixels','Position',[100 100 960 540],'PaperUnits','inches', 'PaperSize',SlideSize, 'PaperPosition',[0 0 SlideSize], ...
    'PaperPositionMode','manual');

subplot(2,2,1); grid on; box on; hold on; 
plot(abs(ExpFRF.h11(1,:)), abs(ExpFRF.h11(2,:)),'LineWidth',1.5,'DisplayName','g1x1x');
plot(abs(ExpFRF.h12(1,:)), abs(ExpFRF.h12(2,:)),'LineWidth',1.5,'DisplayName','g1x2x');
plot(abs(ExpFRF.h33(1,:)), abs(ExpFRF.h33(2,:)),'LineWidth',1.5,'DisplayName','h3x3x');
title('IRCSA | Experimental FRF with Short and Long Bars'); legend show; 
xlim(FigOpts.Xlim); ylim(FigOpts.Ylim)

subplot(2,2,2); hold on; grid on; box on;
plot(ToolCoupledRCSAFRF.Freq,abs(ToolCoupledRCSAFRF.G11),'Color',[0.93 0.69 0.13],'LineWidth',1.5,'DisplayName','Simulated RCSA')
[~,NPoint] = size(ToolRestrainedFRF.Magnitude);
for idxPoint = 1:NPoint
    plot(ToolRestrainedFRF.Freq,ToolRestrainedFRF.Magnitude(:,idxPoint),'Color',[0.49 0.18 0.56],'LineWidth',1.5,'DisplayName','Simulated (Restrained)')
end
title('RCSA | Tool'); xlabel('Frequency (Hz)'); ylabel('Amplitude (m/N)')
xlim(FigOpts.Xlim);ylim(FigOpts.Ylim) 
legend show
plotbrowser

end % end of function
