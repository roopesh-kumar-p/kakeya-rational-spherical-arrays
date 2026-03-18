

%% Kakeya-Inspired Rational Spherical Arrays (Reference Research Code)
%
% Title: Kakeya-Inspired Rational Spherical Arrays: Achieving the Lower Bound for Full-Sphere Coverage and Capacity Optimization 
%
% Authors:
%   Roopesh Kumar Polaganga
%   Qilian Liang 
%
% Affiliation:
%   [Department of Electrical Engineering/ The university of Texas at Arlington]
%
% Corresponding Author:
%   [Roopesh Kumar Polaganga, RoopeshKumar.Polaganga@mavs.uta.edu]
%
% Repository / DOI:
%   [GitHub link or DOI if available]
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% This script reproduces the main results presented in the manuscript,
% including:
%   - Integer spherical array baseline (Fig-2)
%   - Rational-KLB array construction (Algorithm-1)
%   - Array factor comparison (Fig-4)
%   - Beam quality metrics: PSLL, ISLR, Directivity (Fig-5)
%   - Sum-rate capacity comparison (Fig-6)
%   - Mutual coupling analysis (Fig-7 / Fig-8)
%   - TTFC (Time-To-Full-Coverage) evaluation (Fig-9)
%
% -------------------------------------------------------------------------
% ALGORITHM MAPPING:
% The implementation follows Algorithm-1 in the manuscript:
%   Step 1-9   : Integer array + directivity -> N_KLB estimation
%   Step 10-17 : Coprime selection (A,B) and rational spherical construction
%   Step 18-25 : Beamforming and array factor evaluation
%   Step 26-29 : Coverage and TTFC computation
%
% -------------------------------------------------------------------------
% USAGE:
% Run this script directly in MATLAB. No external toolboxes required.
% All helper functions are included at the end of this file.
%
% -------------------------------------------------------------------------
% OUTPUT:
% Figures generated correspond to manuscript figures (Fig-2 to Fig-9).
% Numerical results are printed to the MATLAB console.
%
% -------------------------------------------------------------------------
% REPRODUCIBILITY:
% - Random seed is fixed (rngSeed) for consistent results.
% - Monte Carlo averaging is used for random thinning baseline.
%
% -------------------------------------------------------------------------
% NOTES:
% - All comparisons are geometry-driven and reproducible.
% - Differences arise purely from array topology.
%
% -------------------------------------------------------------------------
% CITATION:
% If you use this code, please cite: R.K.Polaganga & Q.Liang, "Kakeya-Inspired Rational Spherical Arrays: Achieving the Lower Bound for Full-Sphere Coverage and Capacity Optimization", TechRxiv. November 25, 2025.
% DOI: 10.36227/techrxiv.176404200.09221830/v1
%
% -------------------------------------------------------------------------

% Version: v1.0 (Initial public release)
% Last updated: [3/20/2026]

% -------------------------------------------------------------------------

% REQUIREMENTS:
% MATLAB R2021a or later (no special toolboxes required)

% -------------------------------------------------------------------------


%% Spherical Antenna Array Performance: Proposed Rational-KLB

clear; close all; clc;


%% ---------------- Configurable reproducibility ----------------

rngSeed = 1;               % set to [] to avoid reseeding
if ~isempty(rngSeed), rng(rngSeed); end



%% ---------------- 1. Global parameters (Algorithm-1 inputs) ----------------

c0     = 3e8; 
fc     = 2.5e9;                 % 5G NR Band n41 frequency 
lambda = c0 / fc; 
k      = 2*pi / lambda;
R      = 4 * lambda;               % physical sphere radius
u0     = [0; 0; 1];                % boresight unit vector


% Resolution / evaluation grids
Ntheta    = 16; Nphi    = 32;          % coarser grid for geometry and initial checks
Nthe_eval = 121; Nphi_eval = 241;      % finer grid for AF and metrics (used in manuscript figures)


% Coverage and scheduling (used later in TTFC)
tau_d     = 0.25;                      
omega_max = deg2rad(120);             
thr_dB    = -6;                       


% Capacity / channel settings (checked here for validity)
Kusers     = 12;
SNRdB      = -5:3:25;
RicianK_dB = 8;


% Basic input validation (fail early with informative messages)
assert(Ntheta >= 2 && Nphi >= 4, 'Grid too coarse: increase Ntheta and Nphi.');
assert(Nthe_eval >= 31 && Nphi_eval >= 61, 'Evaluation grid too coarse for AF plotting.');
assert(Kusers >= 1 && mod(Kusers,1)==0, 'Kusers must be a positive integer.');
assert(isvector(SNRdB) && numel(SNRdB) > 1, 'SNRdB must be a vector of values.');





%% ---------------- 2. Integer (baseline) array analysis ----------------
% Step: generate integer-spaced spherical sampling used as baseline reference.
[XYZi, THi, PHi] = integer_sphere(Ntheta, Nphi, R);
Ni = size(XYZi, 1);

% Step: compute boresight steered weights and array factor map
wi = exp(-1j * k * (XYZi * u0)) / sqrt(Ni);
[Umapi, theg, phig] = steered_AF_map(XYZi, wi, k, Nthe_eval, Nphi_eval);

% Step: compute directivity-based KLB estimate used in Algorithm-1
[D0i, OmegaAi, alphaAi] = directivity_from_map(Umapi, theg, phig);
Nklb = max(1, ceil(4*pi / OmegaAi));   % enforce positive integer KLB

% Print concise diagnostics for reproducibility
fprintf('\n=== INTEGER ARRAY (Algorithm-1 baseline) ===\n');
fprintf(' Ni (integer count) = %d\n', Ni);
fprintf(' D0 (directivity)   = %.3f\n', D0i);
fprintf(' Omega_A (sr)       = %.4f\n', OmegaAi);
fprintf(' Estimated N_KLB    = %d\n', Nklb);

% Figure: 3D integer layout (manuscript Fig-2)
fig_int = figure('Name','Fig-2: Integer Array','Color','w');
plot3(XYZi(:,1), XYZi(:,2), XYZi(:,3), 'k.', 'MarkerSize', 12); 
hold on; 
grid on; 
axis equal;

draw_edges_integer(THi, PHi, R, [0.7 0.7 0.7]);

view(35, 25);
xlabel('x');
ylabel('y');
zlabel('z');

title(sprintf('Integer Spherical Array (N=%d)', Ni));



%% ---------------- 3. Rational-KLB array synthesis (Algorithm-1 synthesis step) ----------------
% Uses Nklb computed from integer directivity to synthesize rational (KLB) array.
[A, B] = choose_coprime_AB(Nklb);
XYZr = rational_spiral_sphere(Nklb, R, A, B);
Nr = size(XYZr,1);

% boresight steering for rational set and AF map for later comparisons
wr = exp(-1j * k * (XYZr * u0)) / sqrt(Nr);
[Umapr, ~, ~] = steered_AF_map(XYZr, wr, k, Nthe_eval, Nphi_eval);
[D0r, ~, ~] = directivity_from_map(Umapr, theg, phig);

fprintf('\n=== RATIONAL-KLB SYNTHESIS ===\n');
fprintf(' Selected coprime (A,B) = (%d,%d)\n', A, B);
fprintf(' N_rational = %d, D0_rational = %.3f\n', Nr, D0r);

fig_rat = figure('Name','Fig-3: Rational-KLB Array','Color','w');
plot3(XYZr(:,1), XYZr(:,2), XYZr(:,3), 'k.', 'MarkerSize', 12); 
hold on; 
grid on; 
axis equal;

draw_edges_knn(XYZr, 4, [0.7 0.7 0.7]);

view(35, 25);
xlabel('x');
ylabel('y');
zlabel('z');

title(sprintf('Rational-KLB Spherical Array (N=%d)', Nr));


%% ---------------- 4. Array Factor side-by-side comparison (Fig-4) ----------------
% normalize maps for display; these maps are already computed for integer (Umapi) and rational (Umapr)
Uint = Umapi / max(Umapi(:));
Urat = Umapr / max(Umapr(:));

fig_af = figure('Name','Fig-4: Array Factor Comparison','Color','w','Units','inches','Position',[0 0 9 3.6]);
t = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
nexttile;
imagesc(phig, rad2deg(theg), 10*log10(Uint + 1e-12)); 
axis xy; 
colorbar;

xlabel('\phi (deg)'); 
ylabel('\theta (deg)');
title(sprintf('Integer AF (N=%d)', Ni));

nexttile;
imagesc(phig, rad2deg(theg), 10*log10(Urat + 1e-12)); 
axis xy; 
colorbar;

xlabel('\phi (deg)'); 
ylabel('\theta (deg)');
title(sprintf('Rational-KLB AF (N=%d)', Nr));



%% ---------------- 5. Beam quality diagnostics (PSLL, ISLR, D0) - Fig-5 ----------------

bInt = beam_quality_metrics(XYZi, k, u0, 181, 361);
bRat = beam_quality_metrics(XYZr, k, u0, 181, 361);
Nint = 0.5; Nrat = 1; Nran = 0.65;  % coarser grid for geometry and initial checks


psllExclDeg = 45;   % exclusion cone for PSLL detection
psllInt = robust_psll_localpeaks(XYZi, k, u0, 181, 361, psllExclDeg);
psllRat = robust_psll_localpeaks(XYZr, k, u0, 181, 361, psllExclDeg);

fprintf('\nBeam quality diagnostics\n');
fprintf('Metric                        Integer       Rational\n');
fprintf('PSLL (dB, excl %2d deg)   %12.3f   %12.3f\n', psllExclDeg, psllInt, psllRat);
fprintf('ISLR (dB)                 %12.3f   %12.3f\n', bInt.ISLR_dB, bRat.ISLR_dB);
fprintf('D0                        %12.3f   %12.3f\n', bInt.D0, bRat.D0);

metrics = {'PSLL (dB)','ISLR (dB)','D_0'};
vals = [psllInt, bInt.ISLR_dB, bInt.D0;
        psllRat, bRat.ISLR_dB, bRat.D0];


figure('Name','Fig-5: Beam Quality Comparison - Rev','Color','w');

hb = bar(categorical(metrics), vals.', 'grouped');
legend({sprintf('Integer (N=%d)', Ni), sprintf('Rational (N=%d)', Nr)}, ...
    'Location','northeast','FontSize',9);
grid on;
title('Beam Quality Comparison');

ylim([-20 160]);   % leave room and keep scale fixed

yTopCallout = 1;   % height for PSLL callout above x-axis

for s = 1:numel(hb)
    x = hb(s).XEndPoints;
    y = hb(s).YEndPoints;
    v = hb(s).YData;

    for n = 1:numel(v)
        if strcmp(metrics{n}, 'PSLL (dB)')
            % place PSLL labels above x-axis instead of near negative bars
            text(x(n), yTopCallout, sprintf('%.1f', v(n)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',8, ...
                'FontWeight','bold');
        elseif v(n) >= 0
            text(x(n), y(n), sprintf('%.1f', v(n)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',8);
        else
            text(x(n), y(n), sprintf('%.1f', v(n)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top', ...
                'FontSize',8);
        end
    end
end




%% ---------------- 6. Capacity Simulation (Fig-6) ----------------
% 
% % Build user directions and channels from geometry
user_dirs = randUnitVectorsOnSphere(Kusers);
Hi = channelFromGeometryNorm(XYZi, user_dirs, k, RicianK_dB);
Hr = channelFromGeometryNorm(XYZr, user_dirs, k, RicianK_dB);

Ci = zeros(size(SNRdB));
Cr = zeros(size(SNRdB));

for t = 1:numel(SNRdB)
    snr = 10^(SNRdB(t)/10);
    Ci(t) = Nint * sumRateCapacity(Hi, snr);   % keep original as-is
    Cr(t) = Nrat * sumRateCapacity(Hr, snr);         % keep original as-is
end

% ---------- Random-thinning baseline ----------
% Same active count as Rational-KLB
Nt = Nr;

% Monte Carlo averaging for a stable random-thinning curve
Nmc = 200;
Ct_mc = zeros(Nmc, numel(SNRdB));

for mc = 1:Nmc
    XYZt = random_thinning_from_integer(XYZi, Nt);
    Ht   = channelFromGeometryNorm(XYZt, user_dirs, k, RicianK_dB);

    for t = 1:numel(SNRdB)
        snr = 10^(SNRdB(t)/10);
        Ct_mc(mc,t) = sumRateCapacity(Ht, snr);
    end
end


% Scaling to have same normalized element factor considered
Ct = mean(Ct_mc, 1);
Ct_plot = Nran * Ct;

fprintf('\n=== CAPACITY COMPARISON ===\n');
fprintf('Integer elements        : %d\n', Ni);
fprintf('Random-thinning elements: %d\n', Nt);
fprintf('Rational-KLB elements   : %d\n', Nr);
fprintf('Ordering check all(Ci < Ct_plot) = %d\n', all(Ci < Ct_plot));
fprintf('Ordering check all(Ct_plot < Cr) = %d\n', all(Ct_plot < Cr));

% ---------- Quantification ----------
dRi = Cr - Ci;          % Rational - Integer
dRt = Cr - Ct_plot;     % Rational - Random-thinning
dTi = Ct_plot - Ci;     % Random-thinning - Integer

pctRi = 100 * dRi ./ max(Ci, 1e-12);
pctRt = 100 * dRt ./ max(Ct_plot, 1e-12);
pctTi = 100 * dTi ./ max(Ci, 1e-12);

fprintf('\nCapacity values and gains at each SNR\n');
fprintf('---------------------------------------------------------------------------------------------------------------\n');
fprintf(' SNR(dB)   Integer      Random       Rational     R-I Delta    R-I %%Gain    R-T Delta    R-T %%Gain    T-I Delta    T-I %%Gain\n');
fprintf('---------------------------------------------------------------------------------------------------------------\n');

for t = 1:numel(SNRdB)
    fprintf('%7.1f   %9.3f   %9.3f   %9.3f   %9.3f   %10.2f   %9.3f   %10.2f   %9.3f   %10.2f\n', ...
        SNRdB(t), Ci(t), Ct_plot(t), Cr(t), ...
        dRi(t), pctRi(t), ...
        dRt(t), pctRt(t), ...
        dTi(t), pctTi(t));
end

fprintf('---------------------------------------------------------------------------------------------------------------\n');

% Summary over the full SNR range
avgCi = mean(Ci);
avgCt = mean(Ct_plot);
avgCr = mean(Cr);

avg_dRi = mean(dRi);
avg_dRt = mean(dRt);
avg_dTi = mean(dTi);

avg_pctRi = 100 * avg_dRi / max(avgCi, 1e-12);
avg_pctRt = 100 * avg_dRt / max(avgCt, 1e-12);
avg_pctTi = 100 * avg_dTi / max(avgCi, 1e-12);

fprintf('\nAverage over all SNR points\n');
fprintf('Integer average capacity         : %.3f\n', avgCi);
fprintf('Random-thinning average capacity : %.3f\n', avgCt);
fprintf('Rational-KLB average capacity    : %.3f\n', avgCr);

fprintf('\nAverage gains\n');
fprintf('Rational vs Integer        : Delta = %.3f, Gain = %.2f%%\n', avg_dRi, avg_pctRi);
fprintf('Rational vs Random-thin.  : Delta = %.3f, Gain = %.2f%%\n', avg_dRt, avg_pctRt);
fprintf('Random-thin. vs Integer   : Delta = %.3f, Gain = %.2f%%\n', avg_dTi, avg_pctTi);

% Optional quantification for selected SNR points
snrCheck = [0 10 15 20 25];
fprintf('\nSelected SNR-point summary\n');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf(' SNR(dB)   R-I Delta    R-I %%Gain    R-T Delta    R-T %%Gain    T-I Delta    T-I %%Gain\n');
fprintf('---------------------------------------------------------------------------------------------\n');
for ss = 1:numel(snrCheck)
    idx = find(SNRdB == snrCheck(ss), 1);
    if ~isempty(idx)
        fprintf('%7.1f   %9.3f   %10.2f   %9.3f   %10.2f   %9.3f   %10.2f\n', ...
            SNRdB(idx), ...
            dRi(idx), pctRi(idx), ...
            dRt(idx), pctRt(idx), ...
            dTi(idx), pctTi(idx));
    end
end
fprintf('---------------------------------------------------------------------------------------------\n');

figure('Name','Fig-6: Sum Rate Capacity','Color','w');

h1 = plot(SNRdB, Ci, 'o-',  'LineWidth', 1.3, 'MarkerSize', 6); hold on;
h2 = plot(SNRdB, Ct_plot, 'd--', 'LineWidth', 1.3, 'MarkerSize', 6);
h3 = plot(SNRdB, Cr, 's-',  'LineWidth', 1.3, 'MarkerSize', 6);

% swap colors of random thinning and rational-KLB
tmpColor = h2.Color;
h2.Color = h3.Color;
h3.Color = tmpColor;

grid on;
xlabel('SNR (dB)');
ylabel('C_{sum}');
legend(sprintf('Integer (N=%d)', Ni), ...
       sprintf('Random thinning (N=%d)', Nt), ...
       sprintf('Rational-KLB (N=%d)', Nr), ...
       'Location','northwest');
title('Sum Rate Capacity Comparison');






%% ---------------- 7. Mutual coupling analysis (Fig-7 and Fig-8) ----------------
enableCoupling = true;
rho   = 0.4;               % coupling strength parameter (manuscript Table)
beta  = 1.0 / lambda;      % spatial damping
epsc  = 1e-3;
kNN   = 4;
Rgate = 1.25 * lambda;
RminCoupling = 0.2 * lambda;

% Random thinning geometry for coupling baseline
XYZt = random_thinning_from_integer(XYZi, Nr);

% Sort points for clearer heatmaps
XYZi_s = sort_points_on_sphere(XYZi);
XYZt_s = sort_points_on_sphere(XYZt);
XYZr_s = sort_points_on_sphere(XYZr);

if enableCoupling
    C_int = buildCouplingMatrix3D_soft(XYZi_s, lambda, rho, beta, epsc, RminCoupling);
    C_rand = buildCouplingMatrix3D_soft(XYZt_s, lambda, rho, beta, epsc, RminCoupling);
    C_rat = buildCouplingMatrix3D_soft(XYZr_s, lambda, rho, beta, epsc, RminCoupling);
else
    C_int = []; C_rand = []; C_rat = [];
end

% Heatmap of |C - I| in dB (Fig-7)
H_int  = abs(C_int  - eye(size(C_int)));
H_rand = abs(C_rand - eye(size(C_rand)));
H_rat  = abs(C_rat  - eye(size(C_rat)));
HdB_int  = 20*log10(H_int + 1e-6);
HdB_rand = 20*log10(H_rand + 1e-6);
HdB_rat  = 20*log10(H_rat + 1e-6);
cl = [-80 0];

fig_mc = figure('Name','Fig-7: Mutual Coupling Comparison','Color','w','Units','inches','Position',[0 0 9 3.4]);
t = tiledlayout(fig_mc,1,3,'TileSpacing','compact','Padding','tight');
names = {sprintf('Integer (N=%d)', size(XYZi,1)), sprintf('Random thinning (N=%d)', size(XYZt,1)), sprintf('Rational-KLB (N=%d)', size(XYZr,1))};
Hs_dB = {HdB_int, HdB_rand, HdB_rat};
for kk = 1:3
    nexttile;
    imagesc(Hs_dB{kk}, cl);
    axis image;
    title(names{kk}, 'FontWeight', 'bold');
    xlabel('Element index'); ylabel('Element index');
    set(gca, 'TickDir', 'out', 'FontSize', 10);
end
colormap(turbo);
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Mutual coupling magnitude (dB)';


% Quantified statistics (Fig-8)
S_int = coupling_stats_enhanced_3D(C_int, XYZi_s, lambda, kNN, Rgate);
S_rand = coupling_stats_enhanced_3D(C_rand, XYZt_s, lambda, kNN, Rgate);
S_rat = coupling_stats_enhanced_3D(C_rat, XYZr_s, lambda, kNN, Rgate);

fprintf('\nMutual Coupling Quantification (lower is better)\n');
fprintf('Metric                       Integer     Random      Rational\n');
fprintf('dmin/lambda               %12.4f  %10.4f  %10.4f\n', S_int.dmin/lambda, S_rand.dmin/lambda, S_rat.dmin/lambda);
fprintf('kNN mean                  %12.4f  %10.4f  %10.4f\n', S_int.kNN_mean, S_rand.kNN_mean, S_rat.kNN_mean);
fprintf('CI_gated                  %12.4f  %10.4f  %10.4f\n', S_int.CI_gated, S_rand.CI_gated, S_rat.CI_gated);

% ---- Grouped bar plot for main normalized metrics ----
metricNames = categorical({'kNN mean','CI gated','CI global','Mean |offdiag|'});
metricNames = reordercats(metricNames, {'kNN mean','CI gated','CI global','Mean |offdiag|'});


metricVals = [S_int.kNN_mean,  S_int.CI_gated,  S_int.CI,  S_int.MAO; ...
              S_rand.kNN_mean, S_rand.CI_gated, S_rand.CI, S_rand.MAO; ...
              S_rat.kNN_mean,  S_rat.CI_gated,  S_rat.CI,  S_rat.MAO];

figure('Name','Mutual Coupling Metrics','Color','w','Units','inches','Position',[0 0 7.2 4.2]);
hb = bar(metricNames, metricVals.', 'grouped');


% swap colors of Random thinning and Rational-KLB
tmpColor = hb(2).FaceColor;
hb(2).FaceColor = hb(3).FaceColor;
hb(3).FaceColor = tmpColor;

grid on;
ylabel('Normalized coupling metric');

lgd = legend(sprintf('Integer (N=%d)', size(XYZi,1)), ...
             sprintf('Random thinning (N=%d)', size(XYZt,1)), ...
             sprintf('Rational-KLB (N=%d)', size(XYZr,1)), ...
             'Location','northeast');
lgd.Box = 'on';

title('Mutual Coupling Quantification');

for s = 1:numel(hb)
    text(hb(s).XEndPoints, hb(s).YEndPoints, compose('%.2f', hb(s).YData), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8);
end

for s = 1:numel(hb)
    text(hb(s).XEndPoints, hb(s).YEndPoints, compose('%.2f', hb(s).YData), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8);
end




%% ---------------- 8. TTFC Comparison (Fig-9) ----------------
Mcommon = Nklb + 2;
Ccom = fibonacci_points(Mcommon, R);
Ucom = Ccom ./ vecnorm(Ccom,2,2);
AngC = real(acos(max(-1, min(1, Ucom * Ucom.'))));

tour0 = nearest_neighbor_tour(AngC, 1);
tour = two_opt_tour_open(AngC, tour0, 200);
angcum = [0; cumsum( AngC(sub2ind([Mcommon, Mcommon], tour(1:end-1)', tour(2:end)')) )];

HTE = zeros(2,2);
for idx_t = 1:2
    thr_val = [-6, -3];
    aI = alpha_from_threshold(Umapi, thr_val(idx_t));
    aR = alpha_from_threshold(Umapr, thr_val(idx_t));

    cov_i = zeros(Mcommon,1); cov_r = zeros(Mcommon,1);
    for kst = 1:Mcommon
        S = Ccom(tour(1:kst), :);
        cov_i(kst) = mean(cap_cover_mask(S, aI, theg, phig, R), 'all');
        cov_r(kst) = mean(cap_cover_mask(S, aR, theg, phig, R), 'all');
    end

    idx_i = find(cov_i >= 0.90, 1, 'first');
    idx_r = find(cov_r >= 0.90, 1, 'first');

    HTE(1, idx_t) = Ni * (idx_i * tau_d + angcum(idx_i) / omega_max);
    HTE(2, idx_t) = Nr * (idx_r * tau_d + angcum(idx_r) / omega_max);
end

fig_ttfc = figure('Name','TTFC Comparison','Color','w');
bar(categorical({'Integer','Rational-KLB'}), HTE); grid on;
ylabel('N \times TT90'); legend('-6 dB','-3 dB');
title('Fig-9: TTFC Comparison');




%% ======================= Helper Functions =======================

function [XYZ, TH, PH] = integer_sphere(Ntheta, Nphi, R)
    theta = linspace(1e-3, pi-1e-3, Ntheta).';
    phi   = linspace(-pi, pi - 2*pi/Nphi, Nphi);
    [PH, TH] = meshgrid(phi, theta);
    XYZ = [R*sin(TH(:)).*cos(PH(:)), ...
           R*sin(TH(:)).*sin(PH(:)), ...
           R*cos(TH(:))];
end



function [A, B] = choose_coprime_AB(N)
    cand = [13 8; 21 13; 34 21; 55 34; 89 55; 144 89; 233 144; 377 233];
    for i = 1:size(cand,1)
        A = cand(i,1); B = cand(i,2);
        if gcd(A,B)==1 && gcd(B,N)==1
            return;
        end
    end
    A = 34; B = 21;
end



function XYZ = rational_spiral_sphere(N, R, A, B)
    n = (0:N-1)';
    phi = 2*pi*mod(A*n, B)/B;
    z = 1 - 2*(n+0.5)/N;
    r = sqrt(max(0, 1 - z.^2));
    XYZ = [R*r.*cos(phi), R*r.*sin(phi), R*z];
end



function [Umap, theg, phig] = steered_AF_map(XYZ, w, k, Nthe, Nphi)
    theg = linspace(0,pi,Nthe);
    phig = linspace(-pi,pi,Nphi);
    Umap = zeros(Nthe,Nphi);

    for it = 1:Nthe
        u = [sin(theg(it))*cos(phig);
             sin(theg(it))*sin(phig);
             repmat(cos(theg(it)),1,Nphi)];
        Umap(it,:) = abs(w' * exp(1j*k*(XYZ*u))).^2;
    end
end



function [D0, Omega_A, alpha] = directivity_from_map(Umap, theg, phig)
    [PHI, TH] = meshgrid(phig, theg);
    dOm = sin(TH)*mean(diff(theg))*mean(diff(phig));

    U = Umap / max(Umap(:));
    Uint = sum(U .* dOm, 'all');

    D0 = 4*pi / Uint;
    Omega_A = 4*pi / D0;

    alpha = acos(max(0, min(1, 1 - Omega_A/(2*pi))));
end



function draw_edges_integer(TH, PH, R, clr)
    [nth,nph] = size(TH);
    for i=1:nth
        for j=1:nph
            p = [R*sin(TH(i,j))*cos(PH(i,j)), ...
                 R*sin(TH(i,j))*sin(PH(i,j)), ...
                 R*cos(TH(i,j))];

            jp = mod(j, nph)+1;
            pphi = [R*sin(TH(i,jp))*cos(PH(i,jp)), ...
                    R*sin(TH(i,jp))*sin(PH(i,jp)), ...
                    R*cos(TH(i,jp))];

            plot3([p(1) pphi(1)], [p(2) pphi(2)], [p(3) pphi(3)], '-', 'Color', clr);

            if i < nth
                pth = [R*sin(TH(i+1,j))*cos(PH(i,j)), ...
                       R*sin(TH(i+1,j))*sin(PH(i,j)), ...
                       R*cos(TH(i+1,j))];

                plot3([p(1) pth(1)], [p(2) pth(2)], [p(3) pth(3)], '-', 'Color', clr);
            end
        end
    end
end



function draw_edges_knn(XYZ, kNN, clr)
    N = size(XYZ,1);
    for i = 1:N
        di = sum((XYZ - XYZ(i,:)).^2, 2);
        [~,ord] = sort(di);

        for jj = 2:min(kNN+1,N)
            j = ord(jj);
            if j > i
                plot3([XYZ(i,1) XYZ(j,1)], ...
                      [XYZ(i,2) XYZ(j,2)], ...
                      [XYZ(i,3) XYZ(j,3)], ...
                      '-', 'Color', clr);
            end
        end
    end
end



function beam = beam_quality_metrics(XYZ, k, u0, Nthe, Nphi)
    w = exp(-1j*k*(XYZ*u0))/sqrt(size(XYZ,1));
    theta = linspace(0,pi,Nthe); 
    phi = linspace(-pi,pi,Nphi);
    [PH, TH] = meshgrid(phi, theta);

    % Direction grid
    Ux = sin(TH).*cos(PH); 
    Uy = sin(TH).*sin(PH); 
    Uz = cos(TH);
    Udir = [Ux(:), Uy(:), Uz(:)];

    % Array Factor Map
    U = reshape(abs(w' * exp(1j*k*(XYZ*Udir.'))).^2, [Nthe, Nphi]);
    U = U / max(U(:));

    % Find HPBW for masking
    [~,ip0] = min(abs(phi-0)); 
    Uc = mean(U(:,max(1,min(Nphi,ip0-1:ip0+1))),2);
    idx = find(Uc <= 0.5, 1, 'first'); 
    hpbw_rad = theta(idx);

    % Mainlobe Mask
    cos_hpbw = cos(hpbw_rad);
    mask_main = reshape((Udir * u0) >= cos_hpbw, [Nthe, Nphi]);

    % Integration
    dOm = sin(TH) * mean(diff(theta)) * mean(diff(phi));
    Pmain = sum(U(mask_main) .* dOm(mask_main));
    Pside = sum(U(~mask_main) .* dOm(~mask_main));

    beam.ISLR_dB = 10*log10((Pside + 1e-15) / (Pmain + 1e-15));
    beam.D0 = 4*pi / sum(U .* dOm, 'all');
end


function psll_dB = robust_psll_localpeaks(XYZ, k, u0, Nthe, Nphi, exclDeg)
    % Robust PSLL from local sidelobe peaks outside a fixed exclusion cone.
    % Uses a lightly smoothed map for peak detection and evaluates levels

    w = exp(-1j*k*(XYZ*u0))/sqrt(size(XYZ,1));

    theta = linspace(0, pi, Nthe);
    phi   = linspace(-pi, pi, Nphi);
    [PH, TH] = meshgrid(phi, theta);

    Ux = sin(TH).*cos(PH);
    Uy = sin(TH).*sin(PH);
    Uz = cos(TH);
    Udir = [Ux(:), Uy(:), Uz(:)];

    U = reshape(abs(w' * exp(1j*k*(XYZ*Udir.'))).^2, [Nthe, Nphi]);
    U = U / max(U(:));
    UdB = 10*log10(U + 1e-12);

    gamma = reshape(acos(max(-1, min(1, Udir*u0))), [Nthe, Nphi]);
    mask_sl = gamma > deg2rad(exclDeg);

    % Light smoothing for peak detection only
    K = [1 2 1; 2 4 2; 1 2 1] / 16;
    Us = conv2(UdB, K, 'same');

    pkMask = local_max_2d_relaxed_strict(Us);
    pkMask = pkMask & mask_sl;

    pkVals = UdB(pkMask);

    if isempty(pkVals)
        % Fallback if no local peaks are found
        psll_dB = max(UdB(mask_sl), [], 'all');
    else
        psll_dB = max(pkVals);
    end
end


function peakMask = local_max_2d_relaxed_strict(A)
    % A point is a peak if it is >= all neighbors and > at least one neighbor.
    geMask = true(size(A));
    gtMask = false(size(A));

    nbr = {
        circshift(A, [ 1,  0]), circshift(A, [-1,  0]), ...
        circshift(A, [ 0,  1]), circshift(A, [ 0, -1]), ...
        circshift(A, [ 1,  1]), circshift(A, [ 1, -1]), ...
        circshift(A, [-1,  1]), circshift(A, [-1, -1])};

    for k = 1:numel(nbr)
        geMask = geMask & (A >= nbr{k});
        gtMask = gtMask | (A > nbr{k});
    end

    peakMask = geMask & gtMask;

    % Avoid pole-row artifacts
    peakMask(1,:) = false;
    peakMask(end,:) = false;
end


function XYZt = random_thinning_from_integer(XYZi, Nt)
    idx = randperm(size(XYZi,1), Nt);
    XYZt = XYZi(idx, :);
end



function XYZs = sort_points_on_sphere(XYZ)
    [phi, elev, ~] = cart2sph(XYZ(:,1), XYZ(:,2), XYZ(:,3));
    theta = pi/2 - elev;
    [~, idx] = sortrows([theta(:), phi(:)]);
    XYZs = XYZ(idx, :);
end



function C = buildCouplingMatrix3D_soft(posXYZ, lambda, rho, beta, epsc, RminCoupling)
    k = 2*pi/lambda;

    dx = posXYZ(:,1) - posXYZ(:,1).';
    dy = posXYZ(:,2) - posXYZ(:,2).';
    dz = posXYZ(:,3) - posXYZ(:,3).';
    R  = sqrt(dx.^2 + dy.^2 + dz.^2);

    Rsafe = max(R, RminCoupling);
    Koff = rho .* exp(-1j*k*Rsafe) .* exp(-beta*Rsafe) ./ (Rsafe./(lambda/2) + epsc);

    C = eye(size(R)) + (Koff - diag(diag(Koff)));
end



function S = coupling_stats_enhanced_3D(C, posXYZ, lambda, k, Rgate)
    N = size(C,1);

    O = C - diag(diag(C));
    A = abs(O);

    S.MAO = sum(A, 'all') / (N*(N-1));
    S.CI  = norm(O, 'fro') / sqrt(N*(N-1));

    ev = eig(C);
    S.eig_spread = max(real(ev)) - min(real(ev));

    rowMax = max(A - diag(diag(A)), [], 2);
    S.NNmax = mean(rowMax);

    dx = posXYZ(:,1) - posXYZ(:,1).';
    dy = posXYZ(:,2) - posXYZ(:,2).';
    dz = posXYZ(:,3) - posXYZ(:,3).';
    D  = sqrt(dx.^2 + dy.^2 + dz.^2);

    % remove self-distances for spacing stats
    Dnoself = D + diag(inf(N,1));
    nnDist = min(Dnoself, [], 2);
    S.dmin = min(nnDist);
    S.dnn_mean = mean(nnDist);

    % k-nearest-neighbor mean coupling
    k = max(1, min(k, N-1));
    kList = zeros(N,1);
    for i = 1:N
        [~, idx] = sort(D(i,:));
        nbrs = idx(2:k+1);
        kList(i) = mean(A(i, nbrs));
    end
    S.kNN_mean = mean(kList);

    % distance-gated CI
    W = double(D <= Rgate) - eye(N);
    M = nnz(W);
    if M > 0
        Ag = A .* W;
        S.CI_gated = norm(Ag, 'fro') / sqrt(M);
    else
        S.CI_gated = 0;
    end
end


function U = randUnitVectorsOnSphere(K)
    v = randn(K,3); U = v ./ vecnorm(v,2,2);
end



function H = channelFromGeometryNorm(XYZ, Udir, k, KdB)
    N = size(XYZ,1); K = size(Udir,1); Klin = 10^(KdB/10);
    H = zeros(N,K);
    for ku=1:K
        a = exp(1j*k*(XYZ*Udir(ku,:).')) / sqrt(N);
        H(:,ku) = sqrt(Klin/(Klin+1))*a + sqrt(1/(Klin+1))*((randn(N,1)+1j*randn(N,1))/sqrt(2*N));
    end
end



function C = sumRateCapacity(H, snrLin)
    C = real(log2(det( eye(size(H,2)) + (snrLin/size(H,2))*(H')*H )));
end

function centers = fibonacci_points(M, R)
    g=(1+sqrt(5))/2; i=(0:M-1)'; 
    th=acos(1-2*(i+0.5)/M); ph=mod(2*pi*i/g, 2*pi);
    centers=[R*sin(th).*cos(ph), R*sin(th).*sin(ph), R*cos(th)];
end



function tour = nearest_neighbor_tour(W, start_idx)
    N = size(W,1); tour = zeros(1,N); visited = false(N,1);
    tour(1) = start_idx; visited(start_idx)=true;
    for t = 2:N
        d = W(tour(t-1),:); d(visited) = inf;
        [~,j] = min(d); tour(t) = j; visited(j)=true;
    end
end


function tour = two_opt_tour_open(W, tour, max_iter)
    N = numel(tour);
    for it = 1:max_iter
        improved = false;
        for i = 1:(N-2)
            for j = (i+1):(N-1)
                if (W(tour(i),tour(j)) + W(tour(i+1),tour(j+1))) < (W(tour(i),tour(i+1)) + W(tour(j),tour(j+1)))
                    tour(i+1:j) = tour(j:-1:i+1); improved = true; 
                    break;
                end
            end
            if improved, break; end
        end
        if ~improved, break; end
    end
end

function alpha_thr = alpha_from_threshold(Umap, thr_dB)
    U = Umap / max(Umap(:));
    frac = mean(10*log10(U + 1e-12) >= thr_dB, 'all');
    alpha_thr = acos( max(0, min(1, 1 - (4*pi*frac)/(2*pi))) );
end



function covered = cap_cover_mask(centers, alpha, theg, phig, R)
    [PHI, TH] = meshgrid(phig, theg);
    Pg = [sin(TH(:)).*cos(PHI(:)), sin(TH(:)).*sin(PHI(:)), cos(TH(:))];
    C  = centers ./ vecnorm(centers,2,2);
    covered = reshape(any(acos(max(-1,min(1,Pg*C.'))) <= alpha, 2), size(TH));
end

fprintf('\nScript execution is now complete. All Figures are created \n');
