clear

%Part 1 Create Phase FIG
figure1 = figure;
axes1 = axes('Parent', figure1);
jiange = 400;  % Grid resolution

% Define coordinate ranges
x = -2:4/(jiange-1):2;
y = -2:4/(jiange-1):2;

% Initialize phase arrays
tp = zeros(jiange, jiange);
tp0 = zeros(jiange, jiange);
tpp = zeros(jiange, jiange);

% Calculate phase patterns
for i = 1:jiange
    for j = 1:jiange
        tp(i,j) = sign(cos((x(1,i)+y(1,j))*pi/4) * cos((x(1,i)-y(1,j))*pi/4));
        tp0(i,j) = sign(sin((x(1,i)+y(1,j))*pi/4) * sin((x(1,i)-y(1,j))*pi/4));
    end
end

% Combine patterns to create phase regions
for i = 1:jiange
    for j = 1:jiange
        if tp(i,j) == 1 && tp0(i,j) == 1
            tpp(i,j) = 1;
        elseif tp(i,j) == 1 && tp0(i,j) == -1
            tpp(i,j) = 2;
        elseif tp(i,j) == -1 && tp0(i,j) == 1
            tpp(i,j) = 3;
        elseif tp(i,j) == -1 && tp0(i,j) == -1
            tpp(i,j) = 4;
        end
    end
end

% Define custom colormap
cmap = [0.5, 0.5, 0.5];  % Gray color
custom_colormap = [194/256, 206/256, 220/256;   % Light blue
                   145/256, 173/256, 158/256;   % Light green
                   216/256, 156/256, 122/256;   % Orange
                   0.5, 0, 0.5];               % Purple

% Plot phase surface
colormap(custom_colormap);
caxis([1, 4]);
hold on
surf(x, y, tpp, 'FaceColor', cmap);
xlabel('\delta/\Omega');
ylabel('t_B/\Omega_0');
zlabel('T');
shading interp;  % Remove grid lines, smooth colors

%Part 2: Parameter Setting and Initialization
l = 10;
hen = 16;        % Grid dimension
m = hen * hen;   % Number of samples
A = zeros(2, m); % Coupling coefficients
N = 40;          % Number of momentum points
omg = 20;        % Number of quasi-energy levels
T = 2*pi/4.4;    % Driving period
neng = 2*pi/T;   % Frequency
g = 0.5*pi/T;    % Coupling strength
fangcha = 0.03;  % Variance parameter

% Set up parameter grid
for i = 1:hen
    for j = 1:hen
        A(1, (i-1)*hen+j) = 4*pi*((j)/(hen+1)-0.5)/T;
        A(2, (i-1)*hen+j) = 4*pi*((i)/(hen+1.5)-0.5)/T;
    end
end

% Preallocate arrays
bohanshu1 = zeros(N, 2*(2*omg+1), m);
bohanshu2 = zeros(N, 2*(2*omg+1), m);
bohanshut1 = zeros(N, N, m);
bohanshut2 = zeros(N, N, m);
bohanshut11 = zeros(N, N, m);
bohanshut22 = zeros(N, N, m);
Heff = zeros(2*(2*omg+1), 2*(2*omg+1));
touying0 = zeros(2, 2, N, N, m);

% Pauli matrices
paolix = [0, 1; 1, 0];
paoliy = [0, -1i; 1i, 0];
paoliz = [1, 0; 0, -1];

% Main calculation loop
for i = 1:m
    for k = 1:N
        k0 = k*2*pi/N;
        
        % Define time-dependent Hamiltonians
        H0 = 1*(-A(1,i)*((sin(k0)*paolix)-cos(k0)*paoliy)+2*g*sin(k0)*paoliz+A(2,i)*paoliy)/2;
        H1 = 1*(-A(1,i)*((sin(k0)*paolix)-cos(k0)*paoliy)-A(2,i)*paoliy)/1i/pi;
        
        Hfu1 = -H1;
        H3 = H1/3;
        Hfu3 = -H3;
        H5 = H1/5;
        Hfu5 = -H5;
        H7 = H1/7;
        Hfu7 = -H7;
        
        % Construct extended Hamiltonian
        for j = 1:2*omg+1
            Heff(2*j-1,2*j-1) = H0(1,1) + (j)*2*pi/T + 2*omg;
            Heff(2*j,2*j) = H0(2,2) + (j)*2*pi/T + 2*omg;
            Heff(2*j-1,2*j) = H0(1,2);
            Heff(2*j,2*j-1) = H0(2,1);
        end
        
        % First harmonic coupling
        for j = 1:2*omg
            Heff(2*j-1,2*j+1) = Hfu1(1,1);
            Heff(2*j,2*j+2) = Hfu1(2,2);
            Heff(2*j-1,2*j+2) = Hfu1(1,2);
            Heff(2*j,2*j+1) = Hfu1(2,1);
            Heff(2*j+1,2*j-1) = H1(1,1);
            Heff(2*j+2,2*j) = H1(2,2);
            Heff(2*j+2,2*j-1) = H1(2,1);
            Heff(2*j+1,2*j) = H1(1,2);
        end
        
        % Third harmonic coupling
        for j = 1:2*omg-2
            Heff(2*j-1,2*j+5) = Hfu3(1,1);
            Heff(2*j,2*j+6) = Hfu3(2,2);
            Heff(2*j-1,2*j+6) = Hfu3(1,2);
            Heff(2*j,2*j+5) = Hfu3(2,1);
            Heff(2*j+5,2*j-1) = H3(1,1);
            Heff(2*j+6,2*j) = H3(2,2);
            Heff(2*j+6,2*j-1) = H3(2,1);
            Heff(2*j+5,2*j) = H3(1,2);
        end
        
        % Fifth harmonic coupling
        for j = 1:2*omg-4
            Heff(2*j-1,2*j+9) = Hfu5(1,1);
            Heff(2*j,2*j+10) = Hfu5(2,2);
            Heff(2*j-1,2*j+10) = Hfu5(1,2);
            Heff(2*j,2*j+9) = Hfu5(2,1);
            Heff(2*j+9,2*j-1) = H5(1,1);
            Heff(2*j+10,2*j) = H5(2,2);
            Heff(2*j+10,2*j-1) = H5(2,1);
            Heff(2*j+9,2*j) = H5(1,2);
        end
        
        % Seventh harmonic coupling
        for j = 1:2*omg-6
            Heff(2*j-1,2*j+13) = Hfu7(1,1);
            Heff(2*j,2*j+14) = Hfu7(2,2);
            Heff(2*j-1,2*j+14) = Hfu7(1,2);
            Heff(2*j,2*j+13) = Hfu7(2,1);
            Heff(2*j+13,2*j-1) = H7(1,1);
            Heff(2*j+14,2*j) = H7(2,2);
            Heff(2*j+14,2*j-1) = H7(2,1);
            Heff(2*j+13,2*j) = H7(1,2);
        end
        
        % Diagonalize extended Hamiltonian
        [V, E] = eigs(Heff, 2*(2*omg+1));
        e1 = diag(real(E) - (2*omg+(omg+1)*2*pi/T)*eye(2*(2*omg+1)));
        
        for j = 1:2*(2*omg+1)
            bohanshu1(k,j,i) = V(j, 2*omg+1);
            bohanshu2(k,j,i) = V(j, 2*omg+2);
        end
    end
end

% Time evolution via Fourier synthesis
for i = 1:m
    for k = 1:N
        for t = 1:N
            for j = 1:2*omg+1
                phase = exp(-1i*(j-omg-1)*(neng)*t/N*2*pi/neng);
                bohanshut1(k,t,i) = bohanshut1(k,t,i) + phase * bohanshu1(k,2*j-1,i);
                bohanshut2(k,t,i) = bohanshut2(k,t,i) + phase * bohanshu1(k,2*j,i);
                bohanshut11(k,t,i) = bohanshut11(k,t,i) + phase * bohanshu2(k,2*j-1,i);
                bohanshut22(k,t,i) = bohanshut22(k,t,i) + phase * bohanshu2(k,2*j,i);
            end
        end
    end
end

% Calculate projection operators
for i = 1:m
    for kx = 1:N
        for t = 1:N
            psi1 = [bohanshut1(kx,t,i); bohanshut2(kx,t,i)];
            psi2 = [bohanshut11(kx,t,i); bohanshut22(kx,t,i)];
            touying0(:,:,kx,t,i) = 1 * (psi1 * psi1') - 1 * (psi2 * psi2');
        end
    end
end

% Calculate similarity matrix
neihe = ones(m, m);
for j = 1:m
    for o = 1:m
        for kx = 1:N
            for t = 1:N
                det_val = det(touying0(:,:,kx,t,o) + touying0(:,:,kx,t,j));
                neihe(o,j) = (1 - exp(-abs(det_val)^2/(fangcha^2))) * neihe(o,j);
            end
        end
    end
end

% Normalize similarity matrix
z = zeros(1, m);
for i = 1:m
    z(1,i) = sum(neihe(i,:));
end

P = zeros(m, m);
for i = 1:m
    for j = 1:m
        P(i,j) = neihe(i,j) / sqrt(z(1,i)) / sqrt(z(1,j));
    end
end

% Eigenvalue decomposition
[V, E] = eigs(P, m);
e = diag(E);

% Hierarchical clustering
for i = 1:8
    dataa(:,i) = V(:,i);
end
data = dataa;
ZZ = linkage(data, 'single');
TT = cluster(ZZ, 'maxclust', 8);

% Create cluster map
f = zeros(hen, hen);
for i = 1:hen
    for j = 1:hen
        for cluster_num = 1:25
            if TT((i-1)*hen+j,1) == cluster_num
                f(i,j) = cluster_num;
            end
        end
    end
end

% Prepare normalized coordinates
for i = 1:hen*hen
    B(1,i) = A(1,i)/pi*T;
    B(2,i) = A(2,i)/pi*T;
end

% Plot clusters with different colors and markers
for i = 1:hen
    for j = 1:hen
        switch f(i,j)
            case 1
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [0 1 0], 'MarkerSize', 40);
            case 2
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [0 0 1], 'MarkerSize', 40);
            case 3
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [1 0 1], 'MarkerSize', 40);
            case 4
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '^-', 'MarkerFaceColor', [0 1 0], ...
                    'MarkerEdgeColor', [0 1 0], 'MarkerSize', 10);
            case 5
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '^-', 'color', [0 0 1], ...
                    'MarkerFaceColor', [0 0 1], 'MarkerSize', 10);
            case 6
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '^-', 'color', [1 0 0], ...
                    'MarkerFaceColor', [1 0 0], 'MarkerSize', 10);
            case 7
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '.', 'color', [1 0 0], 'MarkerSize', 40);
            case 8
                plot3(B(1,(i-1)*hen+j), B(2,(i-1)*hen+j), 10, '^-', 'color', [1 0 1], ...
                    'MarkerFaceColor', [1 0 1], 'MarkerSize', 10);
        end
        hold on
    end
end

% Axis formatting
hold on
set(gca, 'FontName', 'Times New Roman', 'fontsize', 32)
axis square
set(gca, 'LineWidth', 5)
box on
set(gca, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5], 'Layer', 'top')
set(axes1, 'TickLabelInterpreter', 'latex', 'XTick', [-2 -1 0 1 2], 'YTick', [-2 -1 0 1 2], 'FontSize', 32)
ax = gca;
ax.XAxis.TickLabelColor = 'k';
ax.YAxis.TickLabelColor = 'k';
xlabel('$J_1T/\pi$', 'Interpreter', 'latex', 'Color', 'k')
ylabel('$J_2T/\pi$', 'Interpreter', 'latex', 'Color', 'k')

% PCA analysis
X = V(:, 1:8);
[coeff, score, latent, tsquared, explained, mu] = pca(X, 'Algorithm', 'svd');
figure;
plot(score(:,1), score(:,2), 'r*', 'MarkerSize', 10);

% Save data
save('E:\mnist\D.mat', 'e', 'score', 'tpp', 'B', 'f', 'TT');