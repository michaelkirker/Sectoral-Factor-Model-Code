%% Plot Contributions to core inflation
function plot_contrib_graphs(info,results,core)




figure('Name','Relative contribution from Tradable and Non-tradable factors')

CPI_ID = 1;

dates = 1:size(results.factors,1);


contrib = repmat(results.floads(CPI_ID,:,:),size(results.factors,1),1).*results.factors;



bar1 = bar(dates,[median(contrib(:,:,1:results.numfact(1)),2) median(contrib(:,:,results.numfact(1)+1:end),2)],'group');
legend('Tradable factor (median) contribution','Non-tradable factor (median) contribution','location','Northwest')
axis tight





set(bar1(1),'FaceColor',[1 0 0],'edgecolor',[1 0 0]);
set(bar1(2),'FaceColor',[0 0 1],'edgecolor',[0 0 1]);

set(gca, ...
  'Box'         , 'off'     , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'LineWidth'   , 1         );


trunkdates = core.dates(end-size(core.pcpis,1)+1:end);
q1bin = rem(trunkdates,1)<0.01;
q1indx = find(q1bin==1);

datelabels = floor(trunkdates(q1indx));
    

    
    set(gca,'Xtick',q1indx,...
        'Xticklabel',datelabels,'TickDir', 'out');
    
    set(gca,'box','off','ygrid','on');





%% Save graph as PDF
%set(gcf,'paperType','A4','paperunits','normalized','paperPosition',[0,0,1,1],'paperOrientation','landscape');
%print(gcf,'-dpdf',['results/contrib_2_core_' info.countrycode]);
