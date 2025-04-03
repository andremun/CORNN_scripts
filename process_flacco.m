addpath ./daviolinplot/

if ~exist("cornn_bbob_ela.csv","file")
    if exist('data','var')
        clear data;
    end
    
    datadir = 'E:/FLACCO_data/';
    filelist = struct2cell(dir(datadir));
    filelist = filelist(1,3:end);
    data(1,:) = readtable([datadir filelist{1}], "TreatAsMissing", 'NA');
    
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    inc = 2;
    for ii=2:length(filelist)
        aux = readtable([datadir filelist{ii}], "TreatAsMissing", 'NA');
        if size(aux,1)==1
            data(inc,:) = aux;
            inc = inc+1;
        else
            data(inc:inc+1,:) = aux;
            inc = inc+2;
        end
    end
    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

    idx = contains(data.Properties.VariableNames,'runtime') | ...
          contains(data.Properties.VariableNames,'fun_evals') | ...
          contains(data.Properties.VariableNames,'basic') | ...
          all(ismissing(data),1);
    data(:,idx) = [];
    
    inc = 1;
    for ii=1:5
        data.fcname = replace(data.fcname,['_R' num2str(ii)],'');
    end
    instancelist = unique(data.fcname);
    
    for ii=1:length(instancelist)
        idx = contains(data.fcname,instancelist{ii});
        if ii == 1
            data_avg = varfun(@mean, data(idx,:), "ErrorHandler", @errorFunc);
            data_avg = convertvars(data_avg, 'mean_fcname', 'string');
        else
            data_avg(ii,:) = varfun(@mean, data(idx,:), "ErrorHandler", @errorFunc);
        end
        data_avg.mean_Var1(ii) = ii;
        data_avg.mean_fcname(ii) = instancelist{ii};
        inc = inc+1;
    end
    data_avg = renamevars(data_avg,data_avg.Properties.VariableNames,data.Properties.VariableNames);
    writetable(data_avg,"cornn_bbob_ela.csv");
end

%%
opts = delimitedTextImportOptions("NumVariables", 70);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "fcname", "ela_distr_skewness", "ela_distr_kurtosis", "ela_distr_number_of_peaks", "ela_level_mmce_lda_10", "ela_level_mmce_qda_10", "ela_level_mmce_mda_10", "ela_level_lda_qda_10", "ela_level_lda_mda_10", "ela_level_qda_mda_10", "ela_level_mmce_lda_25", "ela_level_mmce_qda_25", "ela_level_mmce_mda_25", "ela_level_lda_qda_25", "ela_level_lda_mda_25", "ela_level_qda_mda_25", "ela_level_mmce_lda_50", "ela_level_mmce_qda_50", "ela_level_mmce_mda_50", "ela_level_lda_qda_50", "ela_level_lda_mda_50", "ela_level_qda_mda_50", "ela_meta_lin_simple_adj_r2", "ela_meta_lin_simple_intercept", "ela_meta_lin_simple_coef_min", "ela_meta_lin_simple_coef_max", "ela_meta_lin_simple_coef_max_by_min", "ela_meta_lin_w_interact_adj_r2", "ela_meta_quad_simple_adj_r2", "ela_meta_quad_simple_cond", "ela_meta_quad_w_interact_adj_r2", "disp_ratio_mean_02", "disp_ratio_mean_05", "disp_ratio_mean_10", "disp_ratio_mean_25", "disp_ratio_median_02", "disp_ratio_median_05", "disp_ratio_median_10", "disp_ratio_median_25", "disp_diff_mean_02", "disp_diff_mean_05", "disp_diff_mean_10", "disp_diff_mean_25", "disp_diff_median_02", "disp_diff_median_05", "disp_diff_median_10", "disp_diff_median_25", "limo_avg_length_reg", "limo_avg_length_norm", "limo_length_mean", "limo_ratio_mean", "nbc_nn_nb_sd_ratio", "nbc_nn_nb_mean_ratio", "nbc_nn_nb_cor", "nbc_dist_ratio_coeff_var", "nbc_nb_fitness_cor", "pca_expl_var_cov_x", "pca_expl_var_cor_x", "pca_expl_var_cov_init", "pca_expl_var_cor_init", "pca_expl_var_PC1_cov_x", "pca_expl_var_PC1_cor_x", "pca_expl_var_PC1_cov_init", "pca_expl_var_PC1_cor_init", "ic_h_max", "ic_eps_s", "ic_eps_max", "ic_eps_ratio", "ic_m0"];
opts.SelectedVariableNames = ["fcname", "ela_distr_skewness", "ela_distr_kurtosis", "ela_distr_number_of_peaks", "ela_level_mmce_lda_10", "ela_level_mmce_qda_10", "ela_level_mmce_mda_10", "ela_level_lda_qda_10", "ela_level_lda_mda_10", "ela_level_qda_mda_10", "ela_level_mmce_lda_25", "ela_level_mmce_qda_25", "ela_level_mmce_mda_25", "ela_level_lda_qda_25", "ela_level_lda_mda_25", "ela_level_qda_mda_25", "ela_level_mmce_lda_50", "ela_level_mmce_qda_50", "ela_level_mmce_mda_50", "ela_level_lda_qda_50", "ela_level_lda_mda_50", "ela_level_qda_mda_50", "ela_meta_lin_simple_adj_r2", "ela_meta_lin_simple_intercept", "ela_meta_lin_simple_coef_min", "ela_meta_lin_simple_coef_max", "ela_meta_lin_simple_coef_max_by_min", "ela_meta_lin_w_interact_adj_r2", "ela_meta_quad_simple_adj_r2", "ela_meta_quad_simple_cond", "ela_meta_quad_w_interact_adj_r2", "disp_ratio_mean_02", "disp_ratio_mean_05", "disp_ratio_mean_10", "disp_ratio_mean_25", "disp_ratio_median_02", "disp_ratio_median_05", "disp_ratio_median_10", "disp_ratio_median_25", "disp_diff_mean_02", "disp_diff_mean_05", "disp_diff_mean_10", "disp_diff_mean_25", "disp_diff_median_02", "disp_diff_median_05", "disp_diff_median_10", "disp_diff_median_25", "limo_avg_length_reg", "limo_avg_length_norm", "limo_length_mean", "limo_ratio_mean", "nbc_nn_nb_sd_ratio", "nbc_nn_nb_mean_ratio", "nbc_nn_nb_cor", "nbc_dist_ratio_coeff_var", "nbc_nb_fitness_cor", "pca_expl_var_cov_x", "pca_expl_var_cor_x", "pca_expl_var_cov_init", "pca_expl_var_cor_init", "pca_expl_var_PC1_cov_x", "pca_expl_var_PC1_cor_x", "pca_expl_var_PC1_cov_init", "pca_expl_var_PC1_cor_init", "ic_h_max", "ic_eps_s", "ic_eps_max", "ic_eps_ratio", "ic_m0"];
opts.VariableTypes = ["string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "fcname"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "fcname"], "EmptyFieldRule", "auto");
data_avg = readtable("cornn_bbob_ela.csv", opts);

ela = removevars(data_avg,'limo_avg_length_norm');
ela = removevars(ela,{'pca_expl_var_cov_x',...
                      'pca_expl_var_cor_x',...
                      'pca_expl_var_PC1_cov_x',...
                      'pca_expl_var_PC1_cor_x'});
featnames = ela.Properties.VariableNames;

ela = table2array(ela(:,2:end));
[X,mu,sigma] = zscore(ela);
[coeff,score,~,~,explained] = pca(X);
rng('default');
Y = tsne(score(:,1:20));

class = zeros(size(data_avg,1),1);
for ii=1:24
    class(contains(data_avg.fcname,['F' num2str(ii)])) = ii;
end

isrelu = contains(data_avg.fcname,'relu');
istanh = contains(data_avg.fcname,'tanh');
istrain = contains(data_avg.fcname,'training');
istest = contains(data_avg.fcname,'test');
class =  class + 100.*isrelu + 200.*istanh;
cats = categorical(class);
cats = renamecats(cats,{'100','200'},{'ReLU','Tanh'});

families = [1:5; 6:9 NaN; 10:14; 15:19; 20:24];

clf;
gscatter(Y(:,1),Y(:,2),cats,[],'.+^v');
legend('Location','northeastoutside');
xlabel('z_{1}'); ylabel('z_{2}');
axis square; grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',12);
print(gcf,'-dpng',['tsne_ela_bbob+cornn.png']);

for ii=1:5
    clf;
    auxclass = class.*(any(class==families(ii,:),2) | any(class==[100 200],2));
    nclass = length(unique(auxclass))-1;
    clrs = [0.5 0.5 0.5; jet(nclass)];
    auxcats = categorical(auxclass);
    auxcats = renamecats(auxcats,{'0','100','200'},{'BBOB','ReLU','Tanh'});
    gscatter(Y(:,1),Y(:,2),auxcats,clrs,'.+^v');
    legend('Location','northeastoutside');
    xlabel('z_{1}'); ylabel('z_{2}');
    axis square; grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    print(gcf,'-dpng',['tsne_ela_bbob' num2str(ii) '+cornn.png']);
end

isCORNN = contains(data_avg.fcname,'F_');

type = 0.*(isrelu.*istrain) + 1.*(isrelu.*istest) + 2.*(istanh.*istrain) + 3.*(istanh.*istest);
gscatter(Y(isCORNN,1),Y(isCORNN,2),type(isCORNN));
legend('ReLU Train','ReLU Test','Tanh Train','Tanh Test','Location','northeastoutside');
axis square; grid on;
xlabel('z_{1}'); ylabel('z_{2}');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
print(gcf,'-dpng','tsne_ela_cornn_train_test.png');

regnames = data_avg.fcname(isrelu);
idx = cell2mat(strfind(regnames,"_Net_1"));
regnames = replace(extractBefore(regnames,idx),"_"," ");
regnames = replace(regnames,"F ","");
regnames = categorical(regnames);
unique_regnames = categories(regnames);
nfunctions = length(unique_regnames);
clrs = jet(nfunctions);
gscatter(Y(isrelu,1),Y(isrelu,2),regnames,clrs,'.+^v');
lgd = legend('Location','northeastoutside');
lgd.NumColumns = 2;
axis square; grid on;
xlabel('z_{1}'); ylabel('z_{2}');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
print(gcf,'-dpng',['tsne_ela_cornn_fcn.png']);

train_test_dist = zeros(nfunctions,2);
Yaux = Y(isrelu,:);
for ii=1:nfunctions
    train_test_dist(ii,1) = pdist(Yaux(regnames == unique_regnames(ii),:));
end

Yaux = Y(istanh,:);
for ii=1:nfunctions
    train_test_dist(ii,2) = pdist(Yaux(regnames == unique_regnames(ii),:));
end

clf;
daviolinplot(train_test_dist(:),'groups',[zeros(nfunctions,1); ones(nfunctions,1)],...
    'violinwidth',0,...
                             'boxcolors','same','outliers',0,...
                             'box',1,'boxwidth',0.8,'scatter',2,...
                             'scattersize',15,'jitter',1,'scattercolors','same',...
                             'legend',{'ReLU','Tanh'});
ylabel('Distance between train and test');
legend('Location','northeastoutside');
axis square; grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',12);


%%
rng('default');
eva = evalclusters(X', 'kmeans', 'gap', 'KList', 3:size(X,2), ... % minimum of three features
                       'Distance', 'correlation','SearchMethod','firstmaxse');
clust = bsxfun(@eq, eva.OptimalY, 1:eva.OptimalK);
[rho,pval] = corr(X,Y);
rho(isnan(rho) | (pval>0.05) | abs(rho)<0.5) = 0;
strong = [];
for ii=1:2
    [aux,ind] = sort(abs(rho(:,ii).*clust), 'descend');
    ind(aux==0) = NaN;
    strong = [strong ind(1,:)];
end
strong = unique(strong(~isnan(strong)));
cc = 0.*isrelu + isrelu + 2.*istanh;

%%
h1 = figure;
h2 = figure;
for ii=1:length(strong)
    figure(h1);
    clf;
    scatter(Y(:,1),Y(:,2),8,ela(:,strong(ii)),'filled');
    colorbar('EastOutside'); axis square; grid on;
    title(replace(featnames(strong(ii)+1),"_"," "));
    xlabel('z_{1}'); ylabel('z_{2}');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    print(gcf,'-dpng',['tsne_' featnames{strong(ii)+1} '.png']);

    figure(h2)
    subplot(1,4,ii);
    daviolinplot(ela(:,strong(ii)),'groups',cc,'boxcolors','k','outliers',0,...
                                 'box',0,'boxwidth',0.8,'scatter',2,...
                                 'scattersize',15,'jitter',1,'scattercolors','same',...
                                 'xtlabels',{'BBOB','ReLU','Tanh'});
    ylabel(replace(featnames(strong(ii)+1),"_"," "));
    % legend('Location','northeastoutside');
    % axis square; grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    % print(gcf,'-dpng',['violin_' featnames{strong(ii)+1} '.png']);
end
print(gcf,'-dpng',['violin_strong_features.png']);

