%% disclaimer
% The curvature calculation in this script uses a modified version of a 
% script published on the matlab fileexchange by Alireza Dastan, 
% which in tutn is based on a paper from Meyer et al. from 2003
%
% citations:
% Alireza Dastan (2021). Gaussian and mean curvatures calculation on a
% triangulated 3d surface (https://www.mathworks.com/matlabcentral/
% fileexchange/61136-gaussian-and-mean-curvatures-calculation-on-a-triangulated-3d-surface),
% MATLAB Central File Exchange. Retrieved October 6, 2021. 
% 
% Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). 
% Discrete differential-geometry operators for triangulated 2-manifolds. 
% In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
% Available at : http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.24.3427&rep=rep1&type=pdf
%% start
clear; close all;
tic
results = struct;
sp = struct; % struct holding shared properties
sp.m = struct; % membrane properties
sp.Ibar = struct; % Ibar properties
sp.Nbar = struct; % Nbar properties
sp.actin = struct; % actin properties

%% parameters
        % geometry
        L = 10000;
        dL = 100;
        R = 500;
        dR = 100;
        
        % shared properties
        sp.m.c0 = 0;
        sp.m.rigidity = 10;
        sp.m.stretch_modulus = 0;
        sp.m.dh = 1;
        sp.m.D = 3;
        
        sp.Ibar.area = 50;
        sp.Ibar.c0 = -1/20;
        sp.Ibar.mu = 1;
        sp.Ibar.rigidity = 35;
        sp.Ibar.k1 = 0.11;
        sp.Ibar.k2 = 36;
        sp.Ibar.ssc = 1;
        sp.Ibar.mu0 = 0;
        sp.Ibar.Sstart = 0;
        sp.Ibar.start_saturation = 0;
        
        sp.Nbar.area = 50;
        sp.Nbar.c0 = 1/20;
        sp.Nbar.mu = 1;
        sp.Nbar.rigidity = 35;
        sp.Nbar.k1 = 0.11;
        sp.Nbar.k2 = 36;
        sp.Nbar.ssc = 1;
        sp.Nbar.mu0 = 0;
        sp.Nbar.start_saturation = 0;

        sp.actin.kon = 11.6;
        sp.actin.koff = 1.4;
        sp.actin.dh = 2.7;
        sp.actin.conc = 10;
        sp.actin.start_position = -50;
        
        % numerical parameters
        tend = 10;
        plot_on = true;
        save_on = false;
        plot_everyt = 1;
        save_everyt = 2;
        plt.save_plt = false;

%% input parser

%% dependent parameters
sp.Ibar.start_values = 0;
sp.Nbar.start_values = 0;
sp.Ibar.conc = exp(sp.Ibar.mu-sp.Ibar.mu0)*sp.Ibar.ssc;
sp.Nbar.conc = exp(sp.Nbar.mu-sp.Nbar.mu0)*sp.Nbar.ssc;

%% initiate geometry
fprintf('Initiate geometry... '); %report Progress
[verticies, triangles, sp.TR, sp.gridsize] = initiate_geometry(R,L,dR,dL);
fprintf('done\n'); %report Progress

%% initiate simulation
fprintf('Initiate simulation...'); %report Progress
% initial condition
t = 0;
set(verticies,'rho_actin',sp.actin.start_position);
set(verticies,'Ibar_S',sp.Ibar.start_saturation);
set(verticies,'Ibar_values',sp.Ibar.start_values);
set(verticies,'Nbar_S',sp.Nbar.start_saturation);
set(verticies,'Nbar_values',sp.Nbar.start_values);

% memory allocation
N = length(verticies);

% verticies properties
update(verticies);
area0 = get(verticies,'area');
set(verticies,{'area0'},area0);
update_free_energy(verticies,sp);

% initial transition rates
idx = cell2mat(get(verticies,'id'));
[r_up(idx,1),r_down(idx,1)] = membrane_transition_rates(verticies,sp);
[r_IbarOn(idx,1),r_IbarOff(idx,1),r_NbarOn(idx,1),r_NbarOff(idx,1)] = protein_transition_rates(verticies,sp);
Ibar_S = cell2mat(get(verticies,'Ibar_S'));
Nbar_S = cell2mat(get(verticies,'Nbar_S'));
r_actinOn = sp.actin.kon*sp.actin.conc * Ibar_S .* (1-Nbar_S) + sp.actin.koff;
r_actinOff = sp.actin.koff*ones(N,1);

lower_boundary = 1:2 * sp.gridsize(2);
upper_boundary = (N - lower_boundary(end)):N;

if plot_on
    t_plot = plot_everyt;
    plot_state_3DCylinder(verticies,sp,plt);
else
    t_plot = Inf;
end

% if save_on
%     t_save = 0;
%     idx_save = 1;
%     results.x = x;
% else
%     t_save = Inf;
% end
fprintf('done\n'); %report Progress

%% evolve simulation
fprintf('starting simulation...\n');
fprintf('%8.3f sec/%4d sec simulated', t,tend); %report Progress

while t<tend    % transition matrix
    TM = [r_up,r_down,r_IbarOn,r_IbarOff,r_NbarOn,r_NbarOff,r_actinOn,r_actinOff];
    
    % if conditions
    TM([lower_boundary,upper_boundary],:) = 0; % set all rates outside of domain to zero
    free_space = [verticies.rho]-[verticies.rho_actin];
    idx = free_space < sp.m.dh;
    TM(idx,2) = 0; % if space between actin filament tip and membrane < m.dh, set r_down = 0
    idx = free_space < sp.actin.dh;
    TM(idx,7) = 0; % if space between actin filament tip and membrane < actin.dh, set actin.r_on = 0
    
    % Gillespie algorithm
    [row, clm, idx, dt] = Gillespie_algorithm(TM);
    t = t + dt;
    
    %update states
    V = verticies(row);
    if clm == 1 || clm == 2
        % update membrane position
        V.move((idx(1)-idx(2))*sp.m.dh);
        V.attached_triangles.update;
        update(V.nnk1);
        V.nnk1.update_protein_saturation(sp);
        update_free_energy(V.nnk1,sp);
        idx = [V.nnk1.id];
        [r_IbarOn(idx,1),r_IbarOff(idx,1),r_NbarOn(idx,1),r_NbarOff(idx,1)] = protein_transition_rates(V.nnk1,sp);
        idx = [V.nnk2.id];
        [r_up(idx,1),r_down(idx,1)] = membrane_transition_rates(V.nnk2,sp);
        
    elseif clm == 3 || clm == 4
        % update Ibar saturation
        V.Ibar_values = V.Ibar_values + (idx(3)-idx(4));
        V.update_protein_saturation(sp);
        update_free_energy(V,sp);
        [r_IbarOn(V.id,1),r_IbarOff(V.id,1),r_NbarOn(V.id,1),r_NbarOff(V.id,1)] = protein_transition_rates(V,sp);
        idx = [V.nnk1.id];
        [r_up(idx,1),r_down(idx,1)] = membrane_transition_rates(V.nnk1,sp);
        r_actinOn(V.id) = sp.actin.kon*sp.actin.conc * V.Ibar_S .* (1-V.Nbar_S) + sp.actin.koff;
        
    elseif clm == 5 || clm == 6
        % update Nbar saturation
        V.Nbar_values = V.Nbar_values + (idx(5)-idx(6));
        V.update_protein_saturation(sp);
        update_free_energy(V,sp);
        [r_IbarOn(V.id,1),r_IbarOff(V.id,1),r_NbarOn(V.id,1),r_NbarOff(V.id,1)] = protein_transition_rates(V,sp);
        idx = [V.nnk1.id];
        [r_up(idx,1),r_down(idx,1)] = membrane_transition_rates(V.nnk1,sp);
        r_actinOn(V.id) = sp.actin.kon*sp.actin.conc * V.Ibar_S .* (1-V.Nbar_S) + sp.actin.koff;
        
    elseif clm == 7 || clm == 8
        % update actin filament position
        V.rho_actin = V.rho_actin + (idx(7)-idx(8))*sp.actin.dh;
    end
    
    if t>=t_plot
        plot_state_3DCylinder(verticies,sp,plt);
        t_plot = t_plot+plot_everyt;
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8.3f sec/%4d sec simulated', t,tend); %report Progress
end
fprintf('\nsimulation complete\nElapsed time %4f sec\n',round(toc,7)); %report Progress

%% support functions
function [verticies, triangles, TR, sz] = initiate_geometry(R,L,dR,dL)
    % geometry
    [theta,rho,Z,sz] = cylinder_pol(R,L,dR,dL);
    for i = 1:length(Z)
        verticies(i,1) = Vertex(i,theta(i),rho(i),Z(i));
    end
    
    % triangles
    points = [verticies.x;verticies.y;verticies.z]';
    DT = delaunayTriangulation(points);
    [T_fb] = freeBoundary(DT);
    
    %remove triangles at top and bottom of cylinder
    TR = triangulation(T_fb,points);
    F = faceNormal(TR);
    ID = F(:,3)==1 | F(:,3)==-1;
    t_idx = TR.ConnectivityList(~ID,1:3);

    TR = triangulation(t_idx,points);
    neighbours_id = neighbors(TR);
    attached_triangles_id = vertexAttachments(TR);
%     attached_verticies_id = edges(TR);

    for i = 1:size(t_idx,1)
        triangles(i,1) = Triangle(i,verticies(t_idx(i,1)),verticies(t_idx(i,2)),verticies(t_idx(i,3)));
    end
    
    for i = 1:size(t_idx,1)
        id = neighbours_id(i,:);
        triangles(i,1).neighbours = triangles(id(~isnan(id)));
    end
    
    for i = 1:length(verticies)
        verticies(i).attached_triangles = triangles(attached_triangles_id{i});
    end
    
    for i = 1:length(verticies)
        verticies(i).nnk1 = unique([verticies(i).attached_triangles.verticies]);
    end
    
    for i = 1:length(verticies)
        verticies(i).nnk2 = unique([verticies(i).nnk1.nnk1]);
    end
    
end

function [theta,rho,z,sz] = cylinder_pol(r,l,dr,dl)
    Nr = round(2*pi*r/dr);
    L = 0:dl:l;
    n_L = length(L);
    theta = (0:1/Nr:(1-1/Nr))*2*pi;
    THETA = repmat(theta,n_L,1);
    RHO = ones(n_L,Nr)*r;
    Z = repmat(L',1,Nr);
    [theta,rho,z] = grid2points(THETA,RHO,Z);
    sz = size(THETA);
end

function [theta,rho,z] = grid2points(THETA,RHO,Z)
    total_length = size(THETA,1)*size(THETA,2);
    linear_index = 1:total_length;
    THETA_t = THETA';
    RHO_t = RHO';
    Z_t = Z';
    theta = THETA_t(linear_index');
    rho = RHO_t(linear_index');
    z = Z_t(linear_index');
end

function update(verticies)
    for i =1:length(verticies)
        [verticies(i).area, verticies(i).mc] = verticies(i).update;
    end
end

function [area,mc] = update2(verticies)
    N = length(verticies);
    area = zeros(N,1);
    mc = zeros(N,1);
    for i =1:N
        [area(i), mc(i)] = verticies(i).update;
    end
end

function update_free_energy(verticies,sp)
    N = length(verticies);
    if N>1
        c = cell2mat(get(verticies,'mc'));
        area = cell2mat(get(verticies,'area'));
        area0 = cell2mat(get(verticies,'area0'));
        Ibar_S = cell2mat(get(verticies,'Ibar_S'));
        Nbar_S = cell2mat(get(verticies,'Nbar_S'));
        fe = free_energy(c,area,area0,Ibar_S,Nbar_S,sp.Ibar,sp.Nbar,sp.m);
        set(verticies,{'fenergy'},num2cell(fe));
    else
        c = get(verticies,'mc');
        area = get(verticies,'area');
        area0 = get(verticies,'area0');
        Ibar_S = get(verticies,'Ibar_S');
        Nbar_S = get(verticies,'Nbar_S');
        fe = free_energy(c,area,area0,Ibar_S,Nbar_S,sp.Ibar,sp.Nbar,sp.m);
        set(verticies,'fenergy',fe);
    end
end

function fe = free_energy(c,area,area0,Ibar_S,Nbar_S,Ibar,Nbar,m)
    mem_bending = m.rigidity/2*(c-m.c0).^2;
    tension = 0.5*m.stretch_modulus*(area-area0).^2./area0;
    Ibar_bending = Ibar.rigidity/2*Ibar_S.*(c-Ibar.c0).^2;
    Nbar_bending = Nbar.rigidity/2*Nbar_S.*(c-Nbar.c0).^2;
    mixing = 1/Ibar.area*(Ibar_S.*log(Ibar_S+1e-99) + Nbar_S.*log(Nbar_S+1e-99) + ...
        (1-Ibar_S-Nbar_S).*log(1-Ibar_S-Nbar_S+1e-99) - ...
        Ibar.mu.*Ibar_S - Nbar.mu.*Nbar_S);
    fe = (mem_bending  + tension + Ibar_bending + Nbar_bending + mixing).*area;
end

function [r_up,r_down] = membrane_transition_rates(verticies,sp)
    N =length(verticies);
    r_up = zeros(N,1);
    r_down = zeros(N,1);
    for i = 1:N
        verticies(i).move(sp.m.dh)
        verticies(i).attached_triangles.update;
        [area_up, mc_up] = update2(verticies(i).nnk1);

        verticies(i).move(-2*sp.m.dh)
        verticies(i).attached_triangles.update;
        [area_down, mc_down] = update2(verticies(i).nnk1);
        
        verticies(i).move(sp.m.dh);
        verticies(i).attached_triangles.update;
                
        area0 = cell2mat(get(verticies(i).nnk1,'area0'));
        Ibar_S = cell2mat(get(verticies(i).nnk1,'Ibar_S'));
        Nbar_S = cell2mat(get(verticies(i).nnk1,'Nbar_S'));
        E0 = sum(cell2mat(get(verticies(i).nnk1,'fenergy')));
        
        E_up = free_energy(mc_up,area_up,area0,Ibar_S,Nbar_S,sp.Ibar,sp.Nbar,sp.m);
        E_down = free_energy(mc_down,area_down,area0,Ibar_S,Nbar_S,sp.Ibar,sp.Nbar,sp.m);
        
        dE_up = sum(E_up) - E0;
        dE_down = sum(E_down) - E0;
        
        dE_up(round(dE_up,10)==0) = 1e-10; %catch exception: prevent division by zero
        dE_down(round(dE_down,10)==0) = 1e-10;
        
        r_up(i,1) = sp.m.D/sp.m.dh^2*dE_up./(exp(dE_up)-1);
        r_down(i,1) = sp.m.D/sp.m.dh^2*dE_down./(exp(dE_down)-1);
        
    end
end

function [Ibar_r_on,Ibar_r_off,Nbar_r_on,Nbar_r_off] = protein_transition_rates(verticies,sp)   
    N = length(verticies);
    if N>1    
        Ibar_values = cell2mat(get(verticies,'Ibar_values'));
        Nbar_values = cell2mat(get(verticies,'Nbar_values'));
        E0 = cell2mat(get(verticies,'fenergy'));
        mc = cell2mat(get(verticies,'mc'));
        area = cell2mat(get(verticies,'area'));
        area0 = cell2mat(get(verticies,'area0'));
    else
        Ibar_values = get(verticies,'Ibar_values');
        Nbar_values = get(verticies,'Nbar_values');
        E0 = get(verticies,'fenergy');
        mc = get(verticies,'mc');
        area = get(verticies,'area');
        area0 = get(verticies,'area0');
    end
    
    Ibar = sp.Ibar;
    Nbar = sp.Nbar;
    m = sp.m;
    
    % states saturations
    Ibar0 = Ibar_values*Ibar.area./area;
    Nbar0 = Nbar_values*Nbar.area./area;
    Ibar_add = (Ibar_values+1)*Ibar.area./area;
    Nbar_add = (Nbar_values+1)*Nbar.area./area;
    Ibar_sub = (Ibar_values-1)*Ibar.area./area;
    Nbar_sub = (Nbar_values-1)*Nbar.area./area;
    
    % check for over- and underflow
    Ibar_overflow = (Ibar_add+Nbar0)>1;
    Nbar_overflow = (Ibar0+Nbar_add)>1;
    Ibar_underflow = Ibar_sub<0;
    Nbar_underflow = Nbar_sub<0;
    
    % states' energy differences
    dE_Iadd = free_energy(mc,area,area0,Ibar_add,Nbar0,Ibar,Nbar,m)-E0;
    dE_Nadd = free_energy(mc,area,area0,Ibar0,Nbar_add,Ibar,Nbar,m)-E0;
%     dE_Isub = free_energy(c,area,Ibar_sub,Nbar0,Ibar,Nbar,m);
%     dE_Nsub = free_energy(c,area,Ibar0,Nbar_sub,Ibar,Nbar,m);

    % transition rates
    Ibar_r_on = Ibar.k1*exp(-dE_Iadd)*Ibar.conc.*(1-Ibar0-Nbar0);
    Nbar_r_on = Nbar.k1*exp(-dE_Nadd)*Nbar.conc.*(1-Ibar0-Nbar0);
    Ibar_r_off = Ibar.k2*Ibar0;
    Nbar_r_off = Nbar.k2*Nbar0;
    
    % catch exception: over- and underflow
    Ibar_r_on(Ibar_overflow) = 0;
    Nbar_r_on(Nbar_overflow) = 0;
    Ibar_r_off(Ibar_underflow) = 0;
    Nbar_r_off(Nbar_underflow) = 0;

end

function [row, clm, idx, dt] = Gillespie_algorithm(TM)
    % overall transition rate
    TM_sum = sum(TM,2);
    R = sum(TM_sum);
    
    % normalized cummulative transition matrix
    TM_sum_cumsum = cumsum(TM_sum);
    TM_cumsum = cumsum(TM,2);
    TM_cumsum(2:end,:) = TM_cumsum(2:end,:) + TM_sum_cumsum(1:end-1);
    TM_n = TM_cumsum/R;
    
    % time step and update time
    dt = -log(rand)/R;

    % uniform distributed random number
    r = rand;
    
    % find indices corresponding bucket
    [clm,row] = find(TM_n'>r,1);
    idx_mat = zeros(size(TM));
    idx_mat(row,clm) = 1;
    idx = idx_mat(row,:);
end

function plot_state_3DCylinder(verticies,sp,plt)
    % plot current membrane position    
    ax1 = subplot(2,2,[1 3]);
    x = cell2mat(get(verticies,'x'));
    y = cell2mat(get(verticies,'y'));
    z = cell2mat(get(verticies,'z'));
    Ibar_S = cell2mat(get(verticies,'Ibar_S'));
%     Ibar_S = rand(length(Ibar_S),1);
    Nbar_S = cell2mat(get(verticies,'Nbar_S'));
%     Nbar_S = rand(length(Nbar_S),1);

%     col = [Ibar_S,Ibar_S,Nbar_S];
    col1 = [1-Ibar_S,ones(length(Ibar_S),1),ones(length(Ibar_S),1)];
    col2 = [ones(length(Nbar_S),1),1-Nbar_S,ones(length(Nbar_S),1)];

    h1=trisurf(sp.TR.ConnectivityList,x,y,z);
    hold on
    h2=trisurf(sp.TR.ConnectivityList,x,y,z);
    hold off

    axis equal
    title('membrane position & protein saturation')
    xlabel('width [nm]')
    ylabel('height [nm]')
    zlabel('length [nm]')
    cmap = [ones(101,1),[0:.01:1]',ones(101,1)];
    cmap = [cmap;[1:-.01:0]',ones(101,1),ones(101,1)];
%     cmap = [[0:.01:1]',[0:.01:1]',[1:-.01:0]'];
    colormap(ax1,cmap);
    cb = colorbar('Ticks',[.1,.9],...
         'TickLabels',{'Nbar','Ibar'});
    cb.Label.String = 'Protein saturation [%]';
    axis tight
    % Additional bit to control color of each patch individually
    set(gca,'CLim',[0, 1]);
    set(h1,'FaceVertexCData',col1,'CDataMapping','scaled','EdgeColor',[0 0 0],'FaceAlpha',0.5);
    set(h2,'FaceVertexCData',col2,'CDataMapping','scaled','EdgeColor','none','FaceAlpha',0.5);

    
    % plot current energy states
    c = -.1:.001:.1;
    area = mean([verticies.area0]);
    IbarS = 0:.01:1;
    NbarS = 0:.01:1;
    [C_I,IBARS] = meshgrid(c,IbarS);
    [C_N,NBARS] = meshgrid(c,NbarS);
    
    ax2 = subplot(2,2,2);
        E_Ibar = free_energy(C_I,area,area,IBARS,zeros(size(NBARS)),sp.Ibar,sp.Nbar,sp.m);
        surf(C_I,IBARS,E_Ibar)
        hold on
        [E_minLine,idx]=min(E_Ibar,[],2);
        E_min_abs = min(E_Ibar,[],'all');
        [row,clm] = find(E_Ibar == E_min_abs);
        plot3(c(idx),IbarS,E_minLine)
        plot3(c(clm),IbarS(row),E_min_abs,'g*')
        mc = cell2mat(get(verticies,'mc'));
        E0 = cell2mat(get(verticies,'fenergy'));
        plot3(mc(Ibar_S>=0),Ibar_S(Ibar_S>=0),E0(Ibar_S>=0),'ro')
        hold off
        view(-158,35)
        title('free energy f(c, Ibar)')
        xlabel('curvature [1/nm]')
        ylabel('Ibar saturation [%]')
        zlabel('free energy [kB T]')
        legend('energy landscape','minimum line','absolut minimum','data points')
        colormap(ax2,parula);

    ax4 = subplot(2,2,4);
        E_Nbar = free_energy(C_N,area,area,zeros(size(IBARS)),NBARS,sp.Ibar,sp.Nbar,sp.m);
        surf(C_N,NBARS,E_Nbar)
        hold on
        [E_minLine,idx]=min(E_Nbar,[],2);
        E_min_abs = min(E_Nbar,[],'all');
        [row,clm] = find(E_Nbar == E_min_abs);
        plot3(c(idx),NbarS,E_minLine)
        plot3(c(clm),NbarS(row),E_min_abs,'g*')
        plot3(mc(Nbar_S>=0),Nbar_S(Nbar_S>=0),E0(Nbar_S>=0),'ro')
        hold off
        view(-168,43)
        title('free energy f(c, Nbar)')
        xlabel('curvature [1/nm]')
        ylabel('Nbar saturation [%]')
        zlabel('free energy [kB T]')
        legend('energy landscape','minimum line','absolut minimum','data points')
        colormap(ax4,parula);
        
    drawnow;
    if plt.save_plt
        filename = datestr(now,'mmmm-dd-yyyy_HH-MM-SS-FFF');
        filepath = cd;
        saveas(gcf,filepath+"\"+filename,'png');
    end
end