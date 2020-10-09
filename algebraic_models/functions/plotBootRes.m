function [] = plotBootRes(Xd, param_names, idx_fig)

n_param = size(Xd, 1);

num = 1:n_param^2;
idx_mat = reshape(num, n_param, n_param)';
idx_diag = diag(idx_mat);

figure(idx_fig)

for i = 1:n_param
    
    subplot(n_param, n_param, idx_diag(i))
    histogram(Xd(i,:))
    xlabel(param_names{i})
    
    idx_dummy = 1;
    
    for j = i+1:n_param
        subplot(n_param, n_param, idx_diag(i) + idx_dummy)
        scatter(Xd(i,:), Xd(j,:))
        xlabel(param_names{i})
        ylabel(param_names{j})
        
        idx_dummy = idx_dummy + 1;
    end
    
    idx_dummy = 1;
    if (i-1)>0
        for j = i-1:-1:1
            subplot(n_param, n_param, idx_diag(i) - idx_dummy)
            text(0.5, 0.5, num2str(round(corr(Xd(i,:)', Xd(j,:)'),2)));
            axis off
            %plot(1,1,'o','linewidth',2)
            idx_dummy = idx_dummy + 1;
        end
    end
    
end


end

