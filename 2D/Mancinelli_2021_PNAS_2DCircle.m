function Mancinelli_2021_PNAS_2DCircle(varargin)
%% start
tic
m = struct;
Ibar = struct;
Nbar = struct;
actin = struct;
results = struct;

%% biological parameters
m.R = 1000;
m.dR = 15;
m.width = 50;
m.dh = 1;
m.D = 3;
m.stretch_modulus = 0;
m.rigidity = 10;
Ibar.rigidity = 35;
Nbar.rigidity = 35;
m.c0 = 0;
Ibar.c0 = -1/20;
Nbar.c0 = 1/20;
Ibar.k1 = 0.11;
Nbar.k1 = 0.11;
Ibar.k2 = 36;
Nbar.k2 = 36;
% Ibar.conc = 10;
% Nbar.conc = 10;
Ibar.standard_state_concentration = 1;
Nbar.standard_state_concentration = 1;
Ibar.mu0 = 0;
Nbar.mu0 = 0;
Ibar.area = 50;
Nbar.area = 50;
Ibar.mu = 1;
Nbar.mu = 1;
actin.k_on = 11.6;
actin.k_off = 1.4;
% actin.k_on = 0;
% actin.k_off = 0;
actin.conc = 10;
actin.dh = 2.7;
actin.start_position = -50;
Ibar_start_saturation = 0;
Nbar_start_saturation = 0;

%% numerical parameters
tend = 60;
plot_on = true;
save_on = false;
plot_everyt = 2;
save_everyt = 2;
plt.save_plt = false;

%% input parser
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'R'
            m.R = varargin{i+1};
        case 'dR'
            m.dR = varargin{i+1};
        case 'dh'
            m.dh = varargin{i+1};
        case 'D'
            m.D = varargin{i+1};
        case 'membrane stretch modulus'
            m.stretch_modulus = varargin{i+1};
        case 'membrane rigidity'
            m.rigidity = varargin{i+1};
        case 'Ibar rigidity'
            Ibar.rigidity  = varargin{i+1};
        case 'Nbar rigidity'
            Nbar.rigidity = varargin{i+1};
        case 'membrane c0'
            m.c0 = varargin{i+1};
        case 'Ibar c0'
            Ibar.c0 = varargin{i+1};
        case 'Nbar c0'
            Nbar.c0 = varargin{i+1};
        case 'Ibar k1'
            Ibar.k1 = varargin{i+1};
        case 'Nbar k1'
            Nbar.k1 = varargin{i+1};
        case 'Ibar k2'
            Ibar.k2 = varargin{i+1};
        case 'Nbar k2'
            Nbar.k2 = varargin{i+1};
        case 'Ibar standard state concentration'
            Ibar.standard_state_concentration  = varargin{i+1};
        case 'Nbar standard state concentration'
            Nbar.standard_state_concentration = varargin{i+1};              
        case 'Ibar area'
            Ibar.area = varargin{i+1};
        case 'Nbar area'
            Nbar.area = varargin{i+1};
        case 'Ibar mu'
            Ibar.mu  = varargin{i+1};
        case 'Nbar mu'
            Nbar.mu = varargin{i+1};
        case 'Ibar mu0'
            Ibar.mu0  = varargin{i+1};
        case 'Nbar mu0'
            Nbar.mu = varargin{i+1};
        case 'actin k_on'
            actin.k_on = varargin{i+1};
        case 'actin k_off'
            actin.k_off  = varargin{i+1};
        case 'actin concentration'
            actin.conc = varargin{i+1};         
        case 'actin start position'
            actin.start_position  = varargin{i+1};
        case 'Ibar start saturation'
            Ibar_start_saturation = varargin{i+1};
        case 'Nbar start saturation'
            Nbar_start_saturation = varargin{i+1};
        case 'tend'
            tend = varargin{i+1};
        case 'plot on'
            plot_on = varargin{i+1}; 
        case 'save on'
            save_on = varargin{i+1};   
        case 'plot everyt'
            plot_everyt = varargin{i+1}; 
        case 'save everyt'
            save_everyt = varargin{i+1}; 
        case 'save plot'
            plt.save_plt = varargin{i+1};
        case 'membrane width'
            m.width = varargin{i+1};
        otherwise
            error('input parser: Unexpected Input');
    end
end
Ibar.conc = exp(Ibar.mu-Ibar.mu0)*Ibar.standard_state_concentration;
Nbar.conc = exp(Nbar.mu-Nbar.mu0)*Nbar.standard_state_concentration;

%% initiate geometry
fprintf('Initiate geometry... '); %report Progress
[theta,m.h] = initiate_circle(m.R,m.dR);
N = length(theta);

% memory allocation
m.c = zeros(N,1); m.area = zeros(N,1); m.r_up = zeros(N,1);
m.r_down = zeros(N,1);
%initial conditions
if isnumeric(Ibar_start_saturation)
    if Ibar_start_saturation>=0 && Ibar_start_saturation <=1
        Ibar.S = Ibar_start_saturation * ones(N,1);
    else
        error('Unexpeted input: saturation value must be between 0 and 1')
    end
elseif Ibar_start_saturation == "random"
    Ibar.S = rand(N,1);
else
    error('Unexpeted input: saturation value must be numeric and between 0 and 1')
end
if isnumeric(Nbar_start_saturation)
    if Nbar_start_saturation>=0 && Nbar_start_saturation <=1
        Nbar.S = Nbar_start_saturation * ones(N,1);
    else
        error('Unexpeted input: saturation value must be between 0 and 1')
    end
elseif Nbar_start_saturation == "random"
    Nbar.S = rand(N,1);
else
    error('Unexpeted input: saturation value must be numeric and between 0 and 1')
end
actin.h = (m.R+actin.start_position)*ones(N,1);
fprintf('done\n'); %report Progress

%% initiate simulation
fprintf('Initiate simulation...'); %report Progress
% mambrane patch curvature, area and energy
[x,y] = pol2cart([theta(end);theta;theta(1)],[m.h(end);m.h;m.h(1)]);
[m.c,m.area] = update_patch([x,y],m);
m.area0 = m.area;
[x,y] = pol2cart(theta,m.h);
m.E0 = free_energy(m.c,m.area,m.area0,Ibar.S,Nbar.S,Ibar,Nbar,m);
Ibar.values = round(Ibar.S.*m.area/Ibar.area);
Nbar.values = round(Nbar.S.*m.area/Nbar.area);

% initial transition rates
[Ibar.r_on,Ibar.r_off,Nbar.r_on,Nbar.r_off] = protein_transition_rates(m.c,m.area,m.area0,Ibar.values,Nbar.values,m.E0,Ibar,Nbar,m);
[r_up,r_down] = membrane_transition_rates(theta(1:7),m.h(1:7),m.area0(2:6),Ibar.S(2:6),Nbar.S(2:6),m.E0(2:6),Ibar,Nbar,m);
m.r_up(:) = r_up(2);
m.r_down(:) = r_down(2);
actin.r_on = actin.k_on*actin.conc * Ibar.S .* (1-Nbar.S) + actin.k_off;
actin.r_off = actin.k_off*ones(N,1);

if plot_on
    t_plot = 0;
    plot_state_2DCircle(theta,Ibar,Nbar,m,actin,plt);
else
    t_plot = Inf;
end

if save_on
    t_save = 0;
    idx_save = 1;
    results.theta = theta;
else
    t_save = Inf;
end
fprintf('done\n'); %report Progress
t = 0;

%% evolve simulation
fprintf('starting simulation...\n');
fprintf('%8.3f sec/%4d sec simulated', t,tend); %report Progress
while t<tend
    % transition matrix
    TM = [m.r_up,m.r_down,Ibar.r_on,Ibar.r_off,Nbar.r_on,Nbar.r_off,actin.r_on,actin.r_off];
    
    % if conditions
    idx = (m.h - actin.h) < m.dh;
    TM(idx,2) = 0; % if space between actin filament tip and membrane < m.dh, set r_down = 0
    idx = (m.h - actin.h) < actin.dh;
    TM(idx,7) = 0; % if space between actin filament tip and membrane < actin.dh, set actin.r_on = 0
    
    % Gillespie algorithm
    [row, clm, idx, dt] = Gillespie_algorithm(TM);
    t = t + dt;
    
    %update states
    nnk1 = get_neighbours_idx(row,1,N);
    nnk2 = get_neighbours_idx(row,2,N);
    nnk3 = get_neighbours_idx(row,3,N);
    nnk4 = get_neighbours_idx(row,4,N);
    
    % update membrane position
    if clm == 1 || clm == 2
        m.h(row) = m.h(row)+(idx(1)-idx(2))*m.dh;
        [x(row),y(row)] = pol2cart(theta(row),m.h(row));
        [m.c(nnk1),m.area(nnk1)] = update_patch([x(nnk2),y(nnk2)],m);
        
        Ibar.S(nnk1) = Ibar.values(nnk1)*Ibar.area./m.area(nnk1);
        Nbar.S(nnk1) = Nbar.values(nnk1)*Nbar.area./m.area(nnk1);
        if any((Ibar.S+Nbar.S) > 1)
%             if any((Ibar.S+Nbar.S) > 1)
%                 error('Protein overflow occurred')
%             end
            % catch exception: Protein overflow due to area decrease
            k = find((Ibar.S+Nbar.S) > 1);
            for i = 1: length(k)
                if Ibar.S(k(i))>Nbar.S(k(i))
                    Ibar.values(k(i)) = Ibar.values(k(i))-1;
                    Ibar.S(k(i)) = Ibar.values(k(i))*Ibar.area./m.area(k(i));
                else
                    Nbar.values(i) = Nbar.values(i)-1;
                    Nbar.S(k(i)) = Nbar.values(k(i))*Nbar.area./m.area(k(i));
                end
            end
        end
        
        m.E0(nnk1) = free_energy(m.c(nnk1),m.area(nnk1),m.area0(nnk1),Ibar.S(nnk1),Nbar.S(nnk1),Ibar,Nbar,m);
        [Ibar.r_on(nnk1),Ibar.r_off(nnk1),Nbar.r_on(nnk1),Nbar.r_off(nnk1)] = protein_transition_rates(m.c(nnk1),m.area(nnk1),m.area0(nnk1),Ibar.values(nnk1),Nbar.values(nnk1),m.E0(nnk1),Ibar,Nbar,m);
        [m.r_up(nnk2), m.r_down(nnk2)] = membrane_transition_rates(theta(nnk4),m.h(nnk4),m.area0(nnk3),Ibar.S(nnk3),Nbar.S(nnk3),m.E0(nnk3),Ibar,Nbar,m);
        if ~isreal(m.r_up)|| ~isreal(m.r_down)
            keyboard
        end
    % update Ibar saturation
    elseif clm == 3 || clm == 4
        Ibar.values(row) = Ibar.values(row) + (idx(3)-idx(4));
        Ibar.S(row) = Ibar.values(row)*Ibar.area/m.area(row);
        m.E0(row) = free_energy(m.c(row),m.area(row),m.area0(row),Ibar.S(row),Nbar.S(row),Ibar,Nbar,m);
        [Ibar.r_on(row),Ibar.r_off(row),Nbar.r_on(row),Nbar.r_off(row)] = protein_transition_rates(m.c(row),m.area(row),m.area0(row),Ibar.values(row),Nbar.values(row),m.E0(row),Ibar,Nbar,m);
        [m.r_up(nnk1), m.r_down(nnk1)] = membrane_transition_rates(theta(nnk3),m.h(nnk3),m.area0(nnk2),Ibar.S(nnk2),Nbar.S(nnk2),m.E0(nnk2),Ibar,Nbar,m);
        actin.r_on(row) = actin.k_on * actin.conc * Ibar.S(row) * (1-Nbar.S(row)) + actin.k_off;
    
    % update Nbar saturation
    elseif clm == 5 || clm == 6
        Nbar.values(row) = Nbar.values(row) + (idx(5)-idx(6));
        Nbar.S(row) = Nbar.values(row)*Nbar.area/m.area(row);
        m.E0(row) = free_energy(m.c(row),m.area(row),m.area0(row),Ibar.S(row),Nbar.S(row),Ibar,Nbar,m);
        [Ibar.r_on(row),Ibar.r_off(row),Nbar.r_on(row),Nbar.r_off(row)] = protein_transition_rates(m.c(row),m.area(row),m.area0(row),Ibar.values(row),Nbar.values(row),m.E0(row),Ibar,Nbar,m);
        [m.r_up(nnk1), m.r_down(nnk1)] = membrane_transition_rates(theta(nnk3),m.h(nnk3),m.area0(nnk2),Ibar.S(nnk2),Nbar.S(nnk2),m.E0(nnk2),Ibar,Nbar,m);
        actin.r_on(row) = actin.k_on * actin.conc * Ibar.S(row) * (1-Nbar.S(row)) + actin.k_off;
    
    % update actin filament position
    elseif clm == 7 || clm == 8
        actin.h(row) = actin.h(row) + (idx(7)-idx(8))*actin.dh;

    end
    
    % plot current state
    if t>=t_plot
        plot_state_2DCircle(theta,Ibar,Nbar,m, actin,plt);
        t_plot = t_plot+plot_everyt;
    end
    
    % save current state
    if t>=t_save
        results.h(idx_save,:) = m.h;
        results.IbarS(idx_save,:) = Ibar.S;
        results.NbarS(idx_save,:) = Nbar.S;
        results.ActinH(idx_save,:) = actin.h;
        t_save = t_save + save_everyt;
        idx_save = idx_save+1;
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%8.3f sec/%4d sec simulated', t,tend); %report Progress
end
fprintf('\nsimulation complete\nElapsed time %4f sec\n',round(toc,7)); %report Progress

end

%% support functions
function [theta,rho] = initiate_circle(r, dr)
    n_r = round(2*pi*r/dr);
    theta = [(0:1/n_r:(1-1/n_r))*2*pi]';
    rho = r*ones(length(theta),1);
end

function [c,area] = update_patch(xy,m)
    v = diff(xy,1,1)/2;
    v1 = v(1:end-1,:,:);
    v2 = v(2:end,:,:);
    n1 = [-v1(:,2,:),v1(:,1,:)];
    n2 = [-v2(:,2,:),v2(:,1,:)];
    N =  (n1+n2)./vecnorm(n1+n2,2,2);
    dot_p = dot(n1,n2,2)./(vecnorm(n1,2,2).*vecnorm(n2,2,2));
%     dot_p(dot_p>=1)=1;
    alpha = acos(dot_p);
    n_vec = (v2-v1);
    s = sign(dot(N,n_vec,2));
    L = (vecnorm(v1,2,2)+vecnorm(v2,2,2));
    area = L*m.width;
    
    S = vecnorm(v1+v2,2,2);
%     S = vecnorm(v1,2,2)+vecnorm(v2,2,2);
    c = -s.*alpha./S;
end

function fe = free_energy(c,area,area0,Ibar_S,Nbar_S,Ibar,Nbar,m)
    mem_bending = m.rigidity/2*(c-m.c0).^2;
    tension = 0.5*m.stretch_modulus*(area-area0).^2./area0;
    Ibar_bending = Ibar.rigidity/2*Ibar_S.*(c-Ibar.c0).^2;
    Nbar_bending = Nbar.rigidity/2*Nbar_S.*(c-Nbar.c0).^2;
    mixing = 1/50*(Ibar_S.*log(Ibar_S+1e-99) + Nbar_S.*log(Nbar_S+1e-99) + ...
        (1-Ibar_S-Nbar_S).*log(1-Ibar_S-Nbar_S+1e-99) - ...
        Ibar.mu.*Ibar_S - Nbar.mu.*Nbar_S);
    fe = (mem_bending  + tension + Ibar_bending + Nbar_bending + mixing).*area;
end

function [r_up,r_down] = membrane_transition_rates(theta,rho,area0,Ibar_S,Nbar_S,e0,Ibar,Nbar,m)
    tr = [theta,rho];
    num_states = (length(tr)-4)*2;
    TR = zeros(5,2,num_states);
    E0 = zeros(3,1,num_states);
    IBAR_S = zeros(3,1,num_states);
    NBAR_S = zeros(3,1,num_states);
    AREA0 = zeros(3,1,num_states);
    j=1;
    for i = 3:length(tr)-2
        TR(:,:,j:j+1) = repmat(tr(i-2:i+2,:),1,1,2);
        E0(:,1,j:j+1) = repmat(e0(i-2:i),1,1,2);
        IBAR_S(:,1,j:j+1) = repmat(Ibar_S(i-2:i),1,1,2);
        NBAR_S(:,1,j:j+1) = repmat(Nbar_S(i-2:i),1,1,2);
        AREA0(:,1,j:j+1) = repmat(area0(i-2:i),1,1,2);
        j=j+2;
    end
    TR(3,2,1:2:end-1) = TR(3,2,1:2:end-1)+m.dh;
    TR(3,2,2:2:end) = TR(3,2,2:2:end)-m.dh;
    [X,Y] = pol2cart(TR(:,1,:),TR(:,2,:));
    [MC,AREA] = update_patch([X,Y],m);
    E_prime = free_energy(MC,AREA,AREA0,IBAR_S,NBAR_S,Ibar,Nbar,m);
    DE = (E_prime-E0).*[1/sqrt(2),1,1/sqrt(2)]';
%     DE = (E_prime-E0).*[.5,1,.5]';

    dE = reshape(sum(DE,1),num_states,1,1);
%     dE(dE==0) = 1e-10;
    f = dE./(exp(dE)-1);
%     f(abs(f)==Inf) = 1;
    r = m.D/m.dh^2*f;
    r(isnan(r))=0; %catch exception dE==0: (exp(0)-1)=0 -> div0 results in r = NaN
    r_up = r(1:2:end-1);
    r_down = r(2:2:end);
end

function [Ibar_r_on,Ibar_r_off,Nbar_r_on,Nbar_r_off] = protein_transition_rates(c,area,area0,Ibar_values,Nbar_values,E0,Ibar,Nbar,m)
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
    dE_Iadd = free_energy(c,area,area0,Ibar_add,Nbar0,Ibar,Nbar,m)-E0;
    dE_Nadd = free_energy(c,area,area0,Ibar0,Nbar_add,Ibar,Nbar,m)-E0;
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

function [n_bc] = get_neighbours_idx(idx,range,total_length)
    n_ls=(idx-range):(idx+range);
    n_bc=mod(n_ls,total_length);
    n_bc(n_bc==0)=total_length;
end

function plot_state_2DCircle(theta,Ibar,Nbar,m,actin,plt)
    % plot current membrane position
    [x,y] = pol2cart(theta,m.h);
    ax1 = subplot(2,2,1);
        plot(x,y)
        hold on
        cmap = [[1:-.01:0]',[0:.01:1]',zeros(101,1)];
        col = [Nbar.S,Ibar.S,zeros(length(x),1)];
        scatter(x,y,30,col)
        axis equal
        [x_ac,y_ac] = pol2cart(theta,actin.h);
        if actin.k_on > 0
            scatter(x_ac,y_ac,'r*');
        end
        hold off
        colormap(ax1,cmap);
        cb = colorbar('Ticks',[.1,.9],...
         'TickLabels',{'Nbar','Ibar'});
        cb.Label.String = 'Protein saturation [%]';
        axis tight
        title('membrane position & protein saturation')
        xlabel('width [nm]')
        ylabel('hright [nm]')
        
    % plot transition rates
    ax3 = subplot(2,2,3);
        polarplot(theta,m.r_up,theta,m.r_down,theta,Ibar.r_on,theta,Ibar.r_off,theta,Nbar.r_on,theta,Nbar.r_off,theta,actin.r_on, theta,actin.r_off)
        legend('membrane up','membrane down','Ibar on', 'Ibar off', 'Nbar on', 'Nbar off','actin on', 'actin off')
%         set(colorbar,'visible','off')
%         axis tight
        title('transition rates')

     % plot current energy states
    c = -.1:.001:.1;
    area = mean(m.area0);
    IbarS = 0:.01:1;
    NbarS = 0:.01:1;
    [C_I,IBARS] = meshgrid(c,IbarS);
    [C_N,NBARS] = meshgrid(c,NbarS);
    
    ax2 = subplot(2,2,2);
        E_Ibar = free_energy(C_I,area,area,IBARS,zeros(size(NBARS)),Ibar,Nbar,m);
        surf(C_I,IBARS,E_Ibar)
        hold on
        [E_minLine,idx]=min(E_Ibar,[],2);
        E_min_abs = min(E_Ibar,[],'all');
        [row,clm] = find(E_Ibar == E_min_abs);
        plot3(c(idx),IbarS,E_minLine)
        plot3(c(clm),IbarS(row),E_min_abs,'g*')
        plot3(m.c(Ibar.S>=0),Ibar.S(Ibar.S>=0),m.E0(Ibar.S>=0),'ro')
        hold off
        view(-158,35)
        title('free energy f(c, Ibar)')
        xlabel('curvature [1/nm]')
        ylabel('Ibar saturation [%]')
        zlabel('free energy [kB T]')
        legend('energy landscape','minimum line','absolut minimum','data points')
        colormap(ax2,parula);
        
    ax4 = subplot(2,2,4);
        E_Nbar = free_energy(C_N,area,area,zeros(size(IBARS)),NBARS,Ibar,Nbar,m);
        surf(C_N,NBARS,E_Nbar)
        hold on
        [E_minLine,idx]=min(E_Nbar,[],2);
        E_min_abs = min(E_Nbar,[],'all');
        [row,clm] = find(E_Nbar == E_min_abs);
        plot3(c(idx),NbarS,E_minLine)
        plot3(c(clm),NbarS(row),E_min_abs,'g*')
        plot3(m.c(Nbar.S>=0),Nbar.S(Nbar.S>=0),m.E0(Nbar.S>=0),'ro')
        hold off
        view(-168,43)
        title('free energy f(c, Nbar)')
        xlabel('curvature [1/nm]')
        ylabel('Nbar saturation [%]')
        zlabel('free energy [kB T]')
        legend('energy landscape','minimum line','absolut minimum','data points')
        colormap(ax4,parula);

    if plt.save_plt
        filename = datestr(now,'mmmm-dd-yyyy_HH-MM-SS-FFF');
        filepath = cd;
        saveas(gcf,filepath+"\"+filename,'png');
    end
    drawnow
end