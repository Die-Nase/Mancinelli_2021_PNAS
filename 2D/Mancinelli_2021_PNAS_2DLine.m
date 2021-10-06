function results = Mancinelli_2021_PNAS_2DLine(varargin)
%% start
tic
m = struct;
Ibar = struct;
Nbar = struct;
actin = struct;
results = struct;

%%  parameters
% geometry parameters
shape = 'flat';
shape_parameters = [];
m.L = 1000;
m.dL = 15;
m.width = 50;

% biological parameters
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
Ibar.standard_state_concentration = 1;
Nbar.standard_state_concentration = 1;
Ibar.mu0 = 0;
Nbar.mu0 = 0;
Ibar.mu = 1;
Nbar.mu = 1;
Ibar.area = 50;
Nbar.area = 50;
Ibar.mu = 1;
Nbar.mu = 1;
actin.k_on = 0;
actin.k_off = 0;
actin.conc = 10;
actin.dh = 2.7;
actin.start_position = -50;
Ibar_start_saturation = 0;
Nbar_start_saturation = 0;

% numerical parameters
tend = 60;
plot_on = true;
save_on = false;
plot_everyt = 2;
plt.ylimit = [];
plt.save_plt = false;
save_everyt = 2;

%% input parser
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'starting shape'
            shape = varargin{i+1}{1};
            shape_parameters = varargin{i+1}{2};
        case 'L'
            m.L = varargin{i+1};
        case 'dL'
            m.dL = varargin{i+1};
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
            Nbar.mu0 = varargin{i+1};
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
        case 'actin dh'
            actin.dh = varargin{i+1};
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
        case 'plot ylimit'
            plt.ylimit = varargin{i+1};
        case 'membrane width'
            m.width = varargin{i+1};
        otherwise
            error('input parser: Unexpected Input');
    end
end

%% dependent parameters
Ibar.conc = exp(Ibar.mu-Ibar.mu0)*Ibar.standard_state_concentration;
Nbar.conc = exp(Nbar.mu-Nbar.mu0)*Nbar.standard_state_concentration;

%% initiate geometry
fprintf('Initiate geometry... '); %report Progress
[x,m.h] = initiate_line(m.L,m.dL,shape,shape_parameters);
N = length(x);

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
actin.h = actin.start_position*ones(N,1);
fprintf('done\n'); %report Progress

%% initiate simulation
fprintf('Initiate simulation...'); %report Progress
% mambrane patch curvature, area and energy
[m.c(2:end-1),m.area(2:end-1)] = update_patch([x,m.h],m);
m.c(1) = m.c(2); m.c(end) = m.c(end-1);
m.area(1) = m.area(2); m.area(end) = m.area(end-1);
m.area0 = m.area;
m.E0 = free_energy(m.c,m.area,m.area0,Ibar.S,Nbar.S,Ibar,Nbar,m);
Ibar.values = round(Ibar.S.*m.area/Ibar.area);
Nbar.values = round(Nbar.S.*m.area/Nbar.area);

% initial transition rates
[Ibar.r_on,Ibar.r_off,Nbar.r_on,Nbar.r_off] = protein_transition_rates(m.c,m.area,m.area0,Ibar.values,Nbar.values,m.E0,Ibar,Nbar,m);
[m.r_up(3:N-2),m.r_down(3:N-2)] = membrane_transition_rates(x,m.h,m.area0(2:N-1),Ibar.S(2:N-1),Nbar.S(2:N-1),m.E0(2:N-1),Ibar,Nbar,m);
actin.r_on = actin.k_on*actin.conc * Ibar.S .* (1-Nbar.S) + actin.k_off;
actin.r_off = actin.k_off*ones(N,1);

if plot_on
    t_plot = 0;
    plot_state_2DLine(x,Ibar,Nbar,m,actin,plt);
else
    t_plot = Inf;
end

if save_on
    t_save = 0;
    idx_save = 1;
    results.x = x;
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
    TM([1:4,end-3:end],:) = 0; % set all rates outside of domain to zero
    idx = (m.h - actin.h) < m.dh;
    TM(idx,2) = 0; % if space between actin filament tip and membrane < m.dh, set r_down = 0
    idx = (m.h - actin.h) < actin.dh;
    TM(idx,7) = 0; % if space between actin filament tip and membrane < actin.dh, set actin.r_on = 0
    
    % Gillespie algorithm
    [row, clm, idx, dt] = Gillespie_algorithm(TM);
    t = t + dt;
    
    %update states
    nnk1 = [row-1:row+1];
    nnk2 = [row-2:row+2];
    nnk3 = [row-3:row+3];
    nnk4 = [row-4:row+4];
    if any(isnan(TM))
        keyboard
    end
    % update membrane position
    if clm == 1 || clm == 2
        m.h(row) = m.h(row)+(idx(1)-idx(2))*m.dh;
        [m.c(nnk1),m.area(nnk1)] = update_patch([x(nnk2),m.h(nnk2)],m);
        Ibar.S(nnk1) = Ibar.values(nnk1)*Ibar.area./m.area(nnk1);
        Nbar.S(nnk1) = Nbar.values(nnk1)*Nbar.area./m.area(nnk1);
        if any((Ibar.S+Nbar.S) > 1)
%             error('Protein overflow occurred')
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
        [m.r_up(nnk2), m.r_down(nnk2)] = membrane_transition_rates(x(nnk4),m.h(nnk4),m.area0(nnk3),Ibar.S(nnk3),Nbar.S(nnk3),m.E0(nnk3),Ibar,Nbar,m);
    
    % update Ibar saturation
    elseif clm == 3 || clm == 4
        Ibar.values(row) = Ibar.values(row) + (idx(3)-idx(4));
        Ibar.S(row) = Ibar.values(row)*Ibar.area/m.area(row);
        m.E0(row) = free_energy(m.c(row),m.area(row),m.area0(row),Ibar.S(row),Nbar.S(row),Ibar,Nbar,m);
        [Ibar.r_on(row),Ibar.r_off(row),Nbar.r_on(row),Nbar.r_off(row)] = protein_transition_rates(m.c(row),m.area(row),m.area0(row),Ibar.values(row),Nbar.values(row),m.E0(row),Ibar,Nbar,m);
        [m.r_up(nnk1), m.r_down(nnk1)] = membrane_transition_rates(x(nnk3),m.h(nnk3),m.area0(nnk2),Ibar.S(nnk2),Nbar.S(nnk2),m.E0(nnk2),Ibar,Nbar,m);
        actin.r_on(row) = actin.k_on * actin.conc * Ibar.S(row) * (1-Nbar.S(row)) + actin.k_off;
    
    % update Nbar saturation
    elseif clm == 5 || clm == 6
        Nbar.values(row) = Nbar.values(row) + (idx(5)-idx(6));
        Nbar.S(row) = Nbar.values(row)*Nbar.area/m.area(row);
        m.E0(row) = free_energy(m.c(row),m.area(row),m.area0(row),Ibar.S(row),Nbar.S(row),Ibar,Nbar,m);
        [Ibar.r_on(row),Ibar.r_off(row),Nbar.r_on(row),Nbar.r_off(row)] = protein_transition_rates(m.c(row),m.area(row),m.area0(row),Ibar.values(row),Nbar.values(row),m.E0(row),Ibar,Nbar,m);
        [m.r_up(nnk1), m.r_down(nnk1)] = membrane_transition_rates(x(nnk3),m.h(nnk3),m.area0(nnk2),Ibar.S(nnk2),Nbar.S(nnk2),m.E0(nnk2),Ibar,Nbar,m);
        actin.r_on(row) = actin.k_on * actin.conc * Ibar.S(row) * (1-Nbar.S(row)) + actin.k_off;
    
    % update actin filament position
    elseif clm == 7 || clm == 8
        actin.h(row) = actin.h(row) + (idx(7)-idx(8))*actin.dh;

    end
    
    % plot current state
    if t>=t_plot
        plot_state_2DLine(x,Ibar,Nbar,m,actin,plt);
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

%% support funcitons
function [x,y] = initiate_line(L,dL,varargin)
    L = dL * round(L/dL);
    x = [0:dL:L]';
    for i = 1:2:length(varargin)
        switch varargin{i}
        case 'flat'
            y = zeros(length(x),1);
        case 'gaussian'
            a = varargin{i+1}(1);
            b = varargin{i+1}(2);
            c = varargin{i+1}(3);
            y = a*exp(-(x-b).^2./(2*c^2));
        case 'half_circle'
            radius = varargin{i+1};
            if abs(radius) < L/2
                error('radius must be larger then half of the total length')
            end
            theta = asin((x-L/2)/radius);
            y = cos(theta)*radius;
        otherwise
            error('Unexpected Input');
        end
    end
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
    c = s.*alpha./S;
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

function [r_up,r_down] = membrane_transition_rates(x,y,area0,Ibar_S,Nbar_S,e0,Ibar,Nbar,m)
    xy = [x,y];
    num_states = (length(xy)-4)*2;
    XY = zeros(5,2,num_states);
    E0 = zeros(3,1,num_states);
    IBAR_S = zeros(3,1,num_states);
    NBAR_S = zeros(3,1,num_states);
    AREA0 = zeros(3,1,num_states);
    j=1;
    for i = 3:length(xy)-2
        XY(:,:,j:j+1) = repmat(xy(i-2:i+2,:),1,1,2);
        E0(:,1,j:j+1) = repmat(e0(i-2:i),1,1,2);
        IBAR_S(:,1,j:j+1) = repmat(Ibar_S(i-2:i),1,1,2);
        NBAR_S(:,1,j:j+1) = repmat(Nbar_S(i-2:i),1,1,2);
        AREA0(:,1,j:j+1) = repmat(area0(i-2:i),1,1,2);
        j=j+2;
    end
    XY(3,2,1:2:end-1) = XY(3,2,1:2:end-1)+m.dh;
    XY(3,2,2:2:end) = XY(3,2,2:2:end)-m.dh;
    [MC,AREA] = update_patch(XY,m);
    E_prime = free_energy(MC,AREA,AREA0,IBAR_S,NBAR_S,Ibar,Nbar,m);
    DE = (E_prime-E0).*[1/sqrt(2),1,1/sqrt(2)]';
%     DE = (E_prime-E0).*[1,1,1]';

    dE = reshape(sum(DE,1),num_states,1,1);
    dE(round(dE,10)==0) = 1e-10;
    f = dE./(exp(dE)-1);
%     f(abs(f)==Inf) = 1;
    r = m.D/m.dh^2*f;
    if any(abs(r)== Inf)
        idx = abs(r)== Inf;
        r(idx) = 0;
    end
    if any(isnan(r))
        idx = abs(r)== Inf;
        r(idx) = 0;
    end
%     r(isnan(r))=0; %catch exception dE==0: (exp(0)-1)=0 -> div0 results in r = NaN
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

function plot_state_2DLine(x,Ibar,Nbar,m,actin,plt)
    % plot current membrane position
    ax1 = subplot(2,2,1);
        plot(x,m.h,'k')
        hold on
        cmap = [[0:.01:1]',[0:.01:1]',[1:-.01:0]'];
        col = [Ibar.S,Ibar.S,Nbar.S];
        scatter(x,m.h,30,col)
        if actin.k_on > 0
            bar(x,actin.h,'r','BaseValue',actin.start_position)
        end
        hold off
        colormap(ax1,cmap);
        cb = colorbar('Ticks',[.1,.9],...
         'TickLabels',{'Nbar','Ibar'});
        cb.Label.String = 'Protein saturation [%]';
        axis tight
        if ~isempty(plt.ylimit)
            ylim(plt.ylimit)
        end
        title('membrane position & protein saturation')
        xlabel('width [nm]')
        ylabel('hright [nm]')
        
    % plot transition rates
    ax3 = subplot(2,2,3);
        plot(x,m.r_up, x,m.r_down, x,Ibar.r_on, x,Ibar.r_off, x,Nbar.r_on, x,Nbar.r_off, x,actin.r_on, x,actin.r_off)
        legend('membrane up','membrane down','Ibar on', 'Ibar off', 'Nbar on', 'Nbar off','actin on', 'actin off')
        set(colorbar,'visible','off')
        axis tight
        title('transition rates')
        xlabel('width [nm]')
        ylabel('rate [Hz]')
    
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
        
    drawnow;
    if plt.save_plt
        filename = datestr(now,'mmmm-dd-yyyy_HH-MM-SS-FFF');
        filepath = cd;
        saveas(gcf,filepath+"\"+filename,'png');
    end
end