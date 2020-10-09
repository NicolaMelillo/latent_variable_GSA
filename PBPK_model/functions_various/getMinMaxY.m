function [ tgrid, Ymin, Ymax, Y_mean, Y_perc_sup, Y_perc_inf, Y_std, Y_median, Y_CV ] = getMinMaxY( tall, Yall, perc_inf, perc_sup )

% consider to put out of the plot! calculate outside and plot inside!

s_Yall = size(Yall);
t1 = tall{1};
tmax = max(t1);
tgrid = t1; %0:passo:tmax;
l_tgrid = length(tgrid);

% for the first dataset
Y_i = Yall{1};
t_i = tall{1};
n_states = length(Y_i(1,:));
I_i = findclosest(tgrid,t_i);
Ymin = Y_i(I_i,:);
Ymax = Y_i(I_i,:);
Y_prop = zeros(l_tgrid, n_states, s_Yall(1));

Y_i_prop = zeros(l_tgrid,n_states);
Y_prop(:,:,1) = Ymax;

% other dataset
for i=2:s_Yall(1)
    Y_i = Yall{i};
    t_i = tall{i};
    I_i = findclosest(tgrid,t_i);
    for j=1:l_tgrid
        for k=1:n_states
            Y_i_prop(j,k) = Y_i(I_i(j),k);
            if Y_i_prop(j,k)<Ymin(j,k)
                Ymin(j,k)=Y_i_prop(j,k);
            elseif Y_i_prop(j,k)>Ymax(j,k)
                Ymax(j,k)=Y_i_prop(j,k);
            end
        end
    end
    Y_prop(:,:,i) = Y_i_prop;
end

%%% calculate percentile, median, mean

%tgrid = tall(:,1);
ltgrid = length(tgrid);
Y_mean = zeros(ltgrid,n_states);
Y_perc_sup = zeros(ltgrid,n_states);
Y_perc_inf = zeros(ltgrid,n_states);
Y_std = zeros(ltgrid,n_states);
Y_median = zeros(ltgrid,n_states);
Y_CV = zeros(ltgrid,n_states);

for i=1:n_states
    Y_prop_res = reshape_Y_prop(Y_prop(:,i,:));
    [ Y_mean(:,i), Y_perc_sup(:,i), Y_perc_inf(:,i), Y_std(:,i), Y_median(:,i), Y_CV(:,i)  ] = basicstat(Y_prop_res,perc_inf,perc_sup);
end




end


function Y_res = reshape_Y_prop(Y)

    sy = size(Y);
    Y_res = zeros(sy(1),sy(3));
    
    for i=1:sy(1)
        for j=1:sy(3)
            Y_res(i,j) = Y(i,1,j);
        end
    end
    
end






