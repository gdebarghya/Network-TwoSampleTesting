clc
colorseq = {[0 0 1]; [0 0 0]; [1 0 0]; [0 1 0]; [0 1 1]};
styleseq = {'-'; '-'; '-'; '-'; '-'};
markerseq = {'^'; 'x'; 'o'; 's'; 'v'};
%% Experiment 1a
load('results/expt1a_trial1000.mat') 

figure;
h = length(m_all);
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for i=1:h
    for flag = 0:1
        % gather data
        if flag==0
            data = [power0_OpShuff(:,i) power0_FroShuff(:,i) power0_NorTest(:,i)];
        else
            data = [power1_OpShuff(:,i) power1_FroShuff(:,i) power1_NorTest(:,i)];            
        end 
        
        % plot data
        ax = axes();
        pl = plot(n_all,data,'LineWidth',2,'MarkerSize',10);
        for j=1:length(pl)
            pl(j).Color = colorseq{j};
            pl(j).Marker = markerseq{j};
            pl(j).LineStyle = styleseq{j};
        end
  
        if flag==0
            hold on;
            plot(n_all,0.05*ones(1,length(n_all)),':k','LineWidth',1);
            ax.Position(1) = 0.075; 
            ax.YLim = [-0.0025 0.1]; 
            ax.YTick = 0:0.025:0.1;            
        else
            ax.Position(1) = 0.15 + wid;
            ax.YLim = [-0.025 1]; 
            ax.YTick = 0:0.25:1;
        end    
        ax.FontSize = 14; ax.FontWeight = 'bold'; 
        ax.XLim = [n_all(1) n_all(end)];
        ax.XTick = n_all(1):300:n_all(end);
        ax.Position(2:4) = [(1-i/(7*h)-i*hei) wid hei];
    end
end
lgd = legend(pl,{'Boot-Spectral','Boot-Frobenius','Asymp-Normal'});
lgd.FontSize = 16; lgd.FontName = 'Monospaced';
lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_1a','epsc')
close

%% Experiment 1b
load('results/expt1b_trial1000.mat') 

figure;
h = 3;
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for i=1:h
    for flag = 0:1
        % gather data
        if flag==0
            switch i
                case 1
                    data = power0_OpShuff;
                case 2 
                    data = power0_FroShuff;
                case 3
                    data = power0_NorTest;
            end
        else
            switch i
                case 1
                    data = power1_OpShuff;
                case 2
                    data = power1_FroShuff;
                case 3
                    data = power1_NorTest;
            end
        end 
        
        % plot data
        ax = axes();
        pl = semilogx(rho_all,data,'LineWidth',2,'MarkerSize',10);
        for j=1:length(pl)
            pl(j).Color = colorseq{j};
            pl(j).Marker = markerseq{j};
            pl(j).LineStyle = styleseq{j};
        end

        if flag==0
            hold on;
            plot(rho_all,0.05*ones(1,length(rho_all)),':k','LineWidth',1);
            ax.Position(1) = 0.075; 
            ax.YLim = [-0.0025 0.1]; 
            ax.YTick = 0:0.025:0.1;            
        else
            ax.Position(1) = 0.15 + wid;
            ax.YLim = [-0.025 1]; 
            ax.YTick = 0:0.25:1;
        end    
        ax.FontSize = 14; ax.FontWeight = 'bold'; 
        ax.XLim = [rho_all(1) rho_all(end)];
        ax.XTick = rho_all;
        ax.Position(2:4) = [(1-i/(6*h)-i*hei) wid hei];
    end
end
lgd = legend(pl,{'m = 2','m = 4','m = 6','m = 8','m = 10'});
lgd.FontSize = 16; lgd.FontName = 'Monospaced';
lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_1b','epsc')
close

%% Experiment 1c
load('results/expt1c_trial1000.mat') 

figure;
h = 2;
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for i = 1:h
    for flag = 0:1
        % gather data
        if flag==0
            switch i
                case 1
                    data = power0_ChiTest;
                case 2
                    data = power0_NorTest;
            end
        else
            switch i
                case 1
                    data = power1_ChiTest;
                case 2
                    data = power1_NorTest;
            end
        end
        
        % plot data
        ax = axes();
        pl = plot(n_all,data,'LineWidth',2,'MarkerSize',10);
        for j=1:length(pl)
            pl(j).Color = colorseq{j};
            pl(j).Marker = markerseq{j};
            pl(j).LineStyle = styleseq{j};
        end
        
        if flag==0
            hold on;
            plot(n_all,0.05*ones(1,length(n_all)),':k','LineWidth',1);
            ax.Position(1) = 0.075;
            ax.YLim = [-0.0125 0.5];
            ax.YTick = [0 0.05 0.25 0.5];
        else
            ax.Position(1) = 0.15 + wid;
            ax.YLim = [-0.025 1];
            ax.YTick = 0:0.25:1;
        end
        ax.FontSize = 14; ax.FontWeight = 'bold';
        ax.XLim = [n_all(1) n_all(end)];
        ax.XTick = n_all;
        ax.Position(2:4) = [(1-i/(7*h)-i*hei) wid hei];
    end
end

lgd = legend(pl,{'m = 10','m = 20','m = 50','m = 100','m = 200'});
lgd.FontSize = 16; lgd.FontName = 'Monospaced';
lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_1c','epsc')
close

%% Experiment 2a
load('results/expt2a_trial1000.mat') 

figure;
h = length(r_all);
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for i=1:h
    for flag = 0:1
        % gather data
        if flag==0
            data = [power0_ASEBoot(:,i) power0_AdjBoot(:,i) power0_TWTest(:,i)];
        else
            data = [power1_ASEBoot(:,i) power1_AdjBoot(:,i) power1_TWTest(:,i)];
        end 
        
        % plot data
        ax = axes();
        pl = plot(n_all,data,'-ro','LineWidth',2,'MarkerSize',10);
        for j=1:length(pl)
            pl(j).Color = colorseq{j};
            pl(j).Marker = markerseq{j};
            pl(j).LineStyle = styleseq{j};
        end
  
        if flag==0
            hold on;
            plot(n_all,0.05*ones(1,length(n_all)),':k','LineWidth',1);
            ax.Position(1) = 0.075; 
        else
            ax.Position(1) = 0.15 + wid;
        end
        
        ax.FontSize = 14; ax.FontWeight = 'bold';
        ax.XLim = [n_all(1) n_all(end)];
        ax.XTick = n_all(1):300:n_all(end);
        ax.Position(2:4) = [(1-i/(7*h)-i*hei) wid hei];
        ax.YLim = [-0.025 1]; 
        ax.YTick = 0:0.25:1;


    end
end
lgd = legend(pl,{'Boot-ASE','Boot-EPA','Asymp-TW'});
lgd.FontSize = 16; lgd.FontName = 'Monospaced';
lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_2a','epsc')
close


%% Experiment 2b
load('results/expt2b_trial1000.mat') 

figure;
h = 3; 
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for i=1:h
    for flag = 0:1
        % gather data
        if flag==0
            switch i
                case 1
                    data = power0_ASEBoot;
                case 2 
                    data = power0_AdjBoot;
                case 3
                    data = power0_TWTest;
            end
        else
            switch i
                case 1
                    data = power1_ASEBoot;
                case 2 
                    data = power1_AdjBoot;
                case 3
                    data = power1_TWTest;
            end
        end 
        data = data(:,5:-1:1);
        
        % plot data
        ax = axes();
        pl = plot(n_all,data,'LineWidth',2,'MarkerSize',10);
        for j=1:length(pl)
            pl(j).Color = colorseq{j};
            pl(j).Marker = markerseq{j};
            pl(j).LineStyle = styleseq{j};
        end
        if flag==0
            hold on;
            plot(n_all,0.05*ones(1,length(n_all)),':k','LineWidth',1);
            ax.Position(1) = 0.075; 
        else
            ax.Position(1) = 0.15 + wid;
        end
        
        ax.FontSize = 14; ax.FontWeight = 'bold';
        ax.XLim = [n_all(1) n_all(end)];
        ax.XTick = n_all(1):300:n_all(end);
        ax.Position(2:4) = [(1-i/(6*h)-i*hei) wid hei];
        ax.YLim = [-0.025 1]; 
        ax.YTick = 0:0.25:1;

    end
end
lgd = legend(pl,{'\rho = 4','\rho = 2','\rho = 1','\rho = 0.5','\rho = 0.25'});
lgd.FontSize = 16; lgd.FontName = 'Monospaced';
lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_2b','epsc')
close

%% Experiment oregon_1bc
load('results/expt_oregon1_trial100.mat')

figure;
h = 1;
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for flag = 0:1
    % gather data
    if flag==0
        data = mean(-log(pval_p),3);
        x = p_all;
    else
        data = mean(-log(pval_k),3);
        x = k_all;
    end
    
    % plot data
    ax = axes();
    plot(x,data,'-b','LineWidth',2);
    hold on
    %     pl = plot(x,data(:,1),'-b','LineWidth',2);
    %         pl(2).Color = 'r';
    %         pl(2).Marker = 'x';
    %         pl(2).LineWidth = 2;
    %         pl(2).MarkerSize = 10;
    
    plot(x,-log(0.05)*ones(1,length(x)),':k','LineWidth',1);
    ax.YLim = [1 5];
    ax.YTick = 1:5;
    ax.XLim = [x(1) x(end)];
    
    if flag==0
        ax.Position(1) = 0.075;
    else
        ax.Position(1) = 0.15 + wid;
    end
    ax.FontSize = 14; ax.FontWeight = 'bold';
    ax.Position(2:4) = [(1-1/(7*h)-hei) wid hei];
end

% lgd = legend(pl,{'Each network','Average'});
% lgd.FontSize = 16; lgd.FontName = 'Monospaced';
% lgd.Position(1:2) = [(0.175+2*wid) (ax.Position(2)+1/(10*h))];
saveas(gcf,'results/expt_oregon1','epsc')
close

%% Experiment oregon_2
load('results/expt_oregon2_trial100.mat')

data_all = mean(-log(pval_e),3);
x = e_all;
figure;
h = 1;
figwid = 800;
fighei = 300*h;
set(gcf, 'Position', [100 100 figwid fighei]);
wid = 0.3; hei = 2/(3*h);
for flag = 0:1
    % gather data
    if flag==0
        data = data_all(:,1:2:end);
    else
        data = data_all(:,2:2:end);
    end
    
    % plot data
    ax = axes();
    plot(x,data,'-b','LineWidth',2);
    hold on
    
    plot(x,-log(0.05)*ones(1,length(x)),':k','LineWidth',1);
    ax.YLim = [0 8];
    ax.XLim = [x(1) x(end)];
    
    if flag==0
        ax.Position(1) = 0.075;
    else
        ax.Position(1) = 0.15 + wid;
    end
    ax.FontSize = 14; ax.FontWeight = 'bold';
    ax.Position(2:4) = [(1-1/(7*h)-hei) wid hei];
end

saveas(gcf,'results/expt_oregon2','epsc')
close
