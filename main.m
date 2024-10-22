%% Monte Carlo Simulation 
clear;clc;close all;
%% 3.10 a, b
rd_matrix = zeros(10,1);
A_matrix = zeros(10,1);
%for sim_num = 1:1
    %% Set params 
    layer1 = layer(1, [], [], [], [], []); %air layer
    layer2 = layer(1.37, 0, 0.09, 10, 10, 1); %bio layer
    layer3 = layer(1, [], [], [], [], []); %air layer
    layer_list = [layer1, layer2, layer3];% air layer
    n_packets = 15000;
    w_threshold = 5e-4;
    
    m = 10;
    n_r = 200;
    n_z = 200;
    n_a = 200;
    dr = 0.005;
    dz = 0.005;
    da = (0.5*pi)/n_a;
    bio_num = 1;
    out_grid = grids(n_r, n_z, n_a, dr, dz, da);
    
    %% Start propagation in multi layer tissue
    %init figure
    %figure();
    
    for i = 1:n_packets
        %for imaging
        xTrajectory = [];
        yTrajectory = [];
        zTrajectory = [];
        
        %init packet
        p = packet(0, 0, bio_num, layer_list);
        
        xTrajectory(end+1) = p.x;
        yTrajectory(end+1) = p.y;
        zTrajectory(end+1) = p.z;
    
        while p.dead ~= 1
            if p.s == 0
                p = p.step_sample(layer_list);
            end
            p = p.move_me(layer_list);
            if p.s == 0 %absorb scatter
                [p, dw] = p.absorption(layer_list);
                ri = sqrt(p.x^2 + p.y^2);
                if layer_list(p.layer).g ~= 1 %go on scattering
                    p = p.scattering(layer_list);
                    %update A_rz
                    out_grid = out_grid.update_Arz(ri, p.z, dw);
                else
                    out_grid = out_grid.update_A_not_scatter(ri, dw);%NN - should take a look
                end 

            
            else %reflect transmit
                p = p.reflection_transmission(layer_list);
                ri = sqrt(p.x^2 + p.y^2);
                ai = acos(abs(p.uz));
                if p.update == 1
                    out_grid = out_grid.update_Tra(ri, ai, p.w);
                    p.update = 0;
                elseif p.update == 2
                    out_grid = out_grid.update_Rdra(ri, ai, p.w);
                    p.update = 0;
                end
            end
            if p.dead ~= 1 && p.w <= w_threshold
                p = p.termination(m);
            end
            xTrajectory(end+1) = p.x;
            yTrajectory(end+1) = p.y;
            zTrajectory(end+1) = p.z;
        end
        validIndices = (0 <= zTrajectory & zTrajectory <= 0.1) & ...
                       (-0.5 <= yTrajectory & yTrajectory <= 0.5) & ...
                       (-0.5 <= xTrajectory & xTrajectory <= 0.5);
        
         figure(1)
         plot3(xTrajectory(validIndices), yTrajectory(validIndices), zTrajectory(validIndices));
         hold on;
    end
    hold off
    xlabel('X (cm)');
    ylabel('Y (cm)');
    zlabel('Z (cm)');
    title(['Propagation of ', num2str(n_packets), ' photons']);
    grid on;
    textProps = {sprintf('mu_a = %.2d cm^{-1}', layer2.mu_a), ...
                 sprintf('mu_s = %.2d cm^{-1}', layer2.mu_s), ...
                 sprintf('g = %.2f', layer2.g), ...
                 sprintf('d = %.2f', layer2.z1)};
    annotation('textbox', [0.65, 0.15, 0.9, 0.8], 'String', textProps, 'FitBoxToText', 'on', 'BackgroundColor', 'w');

    %% 3.10 a,b
    %normalize Rd_ra
    out_grid = out_grid.normalize_Rdra(dr, da, n_packets, n_r, n_a);
    out_grid = out_grid.normalize_Arz(dr,dz, n_packets, n_r, n_z);
    out_grid = out_grid.normalize_Tra(dr,da, n_packets, n_r, n_a);
    rd_matrix(sim_num) = out_grid.Rd;
    A_matrix(sim_num) = out_grid.A;

end 
rd_mean = mean(rd_matrix);
rd_err = std(rd_matrix);
A_mean = mean(A_matrix);
A_err = std(A_matrix);

%% Plot 3.10c
figure(3)
Fz_1 = out_grid.Az/layer_list(bio_num+1).mu_a;
Fz_137 = out_grid_137.Az/layer_list(bio_num+1).mu_a;
plot(out_grid.z(2:end-1), Fz_1(2:end-1), '-k', out_grid_137.z(2:end-1), Fz_137(2:end-1), '--k', 'LineWidth', 1.5)
hold on
set(gca, 'FontSize',20)
xlabel('z (cm)', 'FontSize',18);
ylabel('Fluence (-)', 'FontSize',18);
axis([0 1 0 10])
set(get(gca,'ylabel'),'rotation',0, 'FontSize',20)
ax = gca;
ax.XTick=[0 0.2 0.4 0.6 0.8 1];
ax.XTickLabel ={'0','0.2','0.4','0.6','0.8','1'};
axis([0 1 0 10]);
legend('n_{rel}=1','n_{rel}=1.37')
ax.YScale = 'log'
title('Comparison of internal Fluence as a function of depth')
p = polyfit(out_grid.z(2:end-1), Fz_1(2:end-1), 1); % p(1) is the slope, p(2) is the intercept
slope = p(1);
intercept = p(2);
str = sprintf('Slope n = 1.37: %.2f', slope);
text(min(out_grid.z(2:end-1)) + 1, max( Fz_1(2:end-1)) - 1, str, 'FontSize', 12, 'Color', 'red');


%% Plot color map

figure(2)
image(out_grid.A_rz)
colorbar
xlabel('r (mm)')
ylabel('z (mm)')
title('Absorption (1/cm^3)')


