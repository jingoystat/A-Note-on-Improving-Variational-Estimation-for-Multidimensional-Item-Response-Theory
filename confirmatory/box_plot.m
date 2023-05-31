idx = ra_bias ~= 0;
if exist('N_converge','var') == 0
   N_converge = Nrep;
end
rbias = [ra_bias(idx); rb_bias(idx); rcorr_bias(idx)];
bias = [a_bias(idx); b_bias(idx); corr_bias(idx)];
names = [1,2,3];
names = categorical([repmat(names(1),N_converge,1); repmat(names(2),N_converge,1); repmat(names(3),N_converge,1)]);
boxchart(names, bias, 'BoxWidth',0.25);
hold on
boxchart(names, rbias, 'WhiskerLineColor',[0.8500 0.3250 0.0980], 'BoxWidth',0.25);
hold on
yline(0);
legend('IS', 'GVEM', 'Location', 'northeast');
xticks(categorical([1,2,3]));
xticklabels({" \alpha", "\it b", "\it c", "\Sigma_{\theta}"});
ax = gca;
ax.FontSize = 16;
ax.FontName = 'times';
ylabel("Bias");
if within == 1
    model = "within";
else
    model = "between";
end
Title = "K="+domain+", N="+person + ", " +model + " item, "+  "correlation is "+r;
title(Title);
set(gca,'LooseInset',get(gca,'TightInset'));
hold off

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
file_name = domain+"_"+person+"_"+model+"_"+ r + "_bias";
saveas(gcf,file_name + ".pdf");


rrmse = [ra_rmse(idx); rb_rmse(idx); rcorr_rmse(idx)];
rmse = [a_rmse(idx); b_rmse(idx); corr_rmse(idx)];
% names = [1,2,3];
% names = categorical([repmat(names(1),Nrep,1); repmat(names(2),Nrep,1); repmat(names(3),Nrep,1)]);
boxchart(names, rmse, 'BoxWidth',0.25);
hold on
boxchart(names, rrmse, 'WhiskerLineColor',[0.8500 0.3250 0.0980], 'BoxWidth',0.25);
legend('IS', 'GVEM', 'Location', 'northeast');
xticks(categorical([1,2,3]));
xticklabels({" \alpha", "\it b", "\it c", "\Sigma_{\theta}"});
ax = gca;
ax.FontSize = 16;
ax.FontName = 'times';
ylabel("RMSE");
if within == 1
    model = "within";
else
    model = "between";
end
Title = "K="+domain+", N="+person + ", " +model + " item, "+  "correlation is "+r;
title(Title);
set(gca,'LooseInset',get(gca,'TightInset'));
hold off

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)])
file_name = domain+"_"+person+"_"+model+"_"+ r + "_rmse";
saveas(gcf,file_name + ".pdf");


% bias = [ra_bias; rb_bias; rcorr_bias; a_bias; b_bias; corr_bias;nan*zeros(Nrep,1); nan*zeros(Nrep,1); nan*zeros(Nrep,1)];
% params = [1,2,3];
% params = categorical([repmat(params(1),Nrep,1); repmat(params(2),Nrep,1); repmat(params(3),Nrep,1);repmat(params(1),Nrep,1); repmat(params(2),Nrep,1); repmat(params(3),Nrep,1);repmat(params(1),Nrep,1); repmat(params(2),Nrep,1); repmat(params(3),Nrep,1)]);
% methods = [1,2];
% methods = categorical([repmat(methods(1), 3*Nrep, 1); repmat(methods(2), 3*Nrep, 1); repmat(methods(2), 3*Nrep, 1)]);
% boxchart(params, bias, 'GroupByColor', methods);