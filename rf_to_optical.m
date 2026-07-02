% =========================================================================
% RF -> Optical Mapping
%   Band-anchored global alignment + interactive line-to-line refinement.
%   Features: editable comb parameters w/ recompute, mode stepper, live
%   mapping-quality readout, mismatch colour-coding, residual plot, export.
% =========================================================================
clear; clc; close all;

%% 1. CONSTANTS, COMB PARAMETERS & THRESHOLDS
c = 2.99792458e8;
P.f_rep_1 = 200.002499650406e6;
P.f_ceo_1 = -40e6;
P.f_rep_2 = 200.002531049574e6;
P.f_ceo_2 = -39e6;
P.map_order = 1;            % mapping polynomial order: 1=linear, 2=quadratic, 3=cubic
                           % (comb dispersion makes the RF->optical map non-linear;
                           %  selectable live in the UI via the "Map order" popup)

snr_threshold   = 8.0;
min_dip_depth   = 1.5;
min_teeth_width = 4;

%% 2. LOAD HITRAN DATABASES
startFolder = '/Users/weihengwang/Documents/MATLAB/DCS/data/'; % CHANGE THIS <---
fprintf('Loading HITRAN continuous curve...\n');
hitranFile = fullfile(startFolder, 'SpectrMol_Ab.coef.txt');
if ~isfile(hitranFile)
    [fName, pName] = uigetfile(fullfile(startFolder, '*.txt'), '1. Select HITRAN continuous curve');
    if isequal(fName, 0), error('Script aborted.'); end
    hitranFile = fullfile(pName, fName);
end
hitranRaw = readmatrix(hitranFile, 'CommentStyle', '#');
hitran_freq_Hz = c .* hitranRaw(:,1) .* 100;
hitran_shape = (size(hitranRaw, 2) >= 3) * hitranRaw(:,3) + (size(hitranRaw, 2) < 3) * hitranRaw(:,2);
hitran_shape = hitran_shape / max(hitran_shape, [], 'omitnan');

fprintf('Loading HITRAN discrete lines...\n');
lineFile = fullfile(startFolder, 'SpectrMol_line.txt');
if ~isfile(lineFile)
    [fName, pName] = uigetfile(fullfile(startFolder, '*.txt'), '2. Select HITRAN lines');
    if isequal(fName, 0), error('Script aborted.'); end
    lineFile = fullfile(pName, fName);
end
fid = fopen(lineFile, 'rt'); rawText = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
rawLines = rawText{1};
raw_hitran_lines_THz = []; raw_line_intensity = [];
for k = 1:length(rawLines)
    str = strtrim(rawLines{k});
    if isempty(str) || startsWith(str, '#'), continue; end
    tokens = sscanf(str, '%f');
    if length(tokens) >= 3
        raw_hitran_lines_THz(end+1) = (c * tokens(2) * 100) / 1e12;
        raw_line_intensity(end+1) = tokens(3);
    end
end
raw_line_intensity = raw_line_intensity / max(raw_line_intensity, [], 'omitnan');
sig_mask = raw_line_intensity > 0.1;                  % only meaningful lines
sig_hitran_THz = raw_hitran_lines_THz(sig_mask);
sig_intensity  = raw_line_intensity(sig_mask);

%% 3. EXTRACT, FILTER & FIT RF DIPS  (independent of the comb parameters)
fprintf('Extracting RF dip profiles...\n');
folderPath = uigetdir(startFolder, '3. Select the FOLDER containing out_*.txt files');
if isequal(folderPath, 0), error('Script aborted.'); end
rawFiles = dir(fullfile(folderPath, 'out_*.txt'));
fileNames = {rawFiles.name};
fileNums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), fileNames);
[~, sortIdx] = sort(fileNums);
rawFiles = rawFiles(sortIdx);

rf_data_cell = {}; rf_raw_cell = {}; fit_curves_cell = {}; valid_filenames = {};
rf_depths = []; rf_centers_Hz = []; snr_array = []; depth_array = []; width_array = [];

for i = 1:length(rawFiles)
    data = readmatrix(fullfile(folderPath, rawFiles(i).name));
    freq = data(:,1); power = data(:,2);

    % --- Comb-grid tooth detection: estimate the comb period from the strong
    % (unabsorbed) teeth, then take one peak per comb cell. Captures absorbed
    % teeth and rejects between-teeth noise maxima.
    [tp0, ti0] = findpeaks(power, 'MinPeakProminence', 3);
    if numel(tp0) < 5, continue; end
    strong = ti0(tp0 >= median(tp0));
    comb_P = median(diff(strong));
    if ~isfinite(comb_P) || comb_P < 2, comb_P = median(diff(ti0)); end
    [teeth_power, teeth_idx] = findpeaks(power, ...
        'MinPeakDistance', max(1, round(0.6*comb_P)), 'MinPeakProminence', 1.5);
    if numel(teeth_power) < min_teeth_width, continue; end
    teeth_freq = freq(teeth_idx);

    % Curved comb-power envelope (robust upper fit) so the Gaussian envelope's
    % slope/curvature over this small RF window doesn't bias the dip extent or
    % centre. Absorption only pushes teeth DOWN, so the envelope tracks the tops.
    [baseline_vec, env_fun] = comb_envelope(teeth_freq, teeth_power);
    dip_profile = max(baseline_vec - teeth_power, 0);
    max_dip = max(dip_profile);

    ds = sort(dip_profile);                            % robust baseline noise
    noise_floor = max(std(ds(1:max(3, round(0.3*numel(ds))))), 0.01);
    snr = max_dip / noise_floor;

    [max_val, maxIdx] = max(dip_profile);              % dip extent
    ext_thr = max(0.15 * max_val, 5 * noise_floor);
    leftIdx = maxIdx;  while leftIdx > 1 && dip_profile(leftIdx) >= ext_thr, leftIdx = leftIdx - 1; end
    rightIdx = maxIdx; while rightIdx < length(dip_profile) && dip_profile(rightIdx) >= ext_thr, rightIdx = rightIdx + 1; end
    v_idx = leftIdx:rightIdx;

    if snr < snr_threshold || max_dip < min_dip_depth || length(v_idx) < min_teeth_width, continue; end

    t_freq_crop = teeth_freq(v_idx); dip_crop = dip_profile(v_idx);
    tp_crop     = teeth_power(v_idx);

    % Bounded Lorentzian fit (all teeth kept) for the sub-tooth centre + curve.
    lor      = @(p,f) p(1).*(p(3)/2).^2 ./ ((f - p(2)).^2 + (p(3)/2).^2);
    tooth_df = median(diff(teeth_freq));
    f_deep   = teeth_freq(maxIdx);
    span_c   = t_freq_crop(end) - t_freq_crop(1);
    c_lo     = max(t_freq_crop(1),   f_deep - 5*tooth_df);
    c_hi     = min(t_freq_crop(end), f_deep + 5*tooth_df);
    lb = [0,       c_lo, max(tooth_df, eps)];
    ub = [max_dip, c_hi, max(span_c, 2*tooth_df)];
    guess = [max_dip, f_deep, span_c/3];
    fit_params = fminsearch(@(p) lor_cost(p, t_freq_crop, dip_crop, lb, ub), guess, optimset('Display','off'));

    rf_raw_cell{end+1}    = [freq, power];
    rf_data_cell{end+1}   = [t_freq_crop, dip_crop, tp_crop, baseline_vec(v_idx)];
    rf_centers_Hz(end+1)  = fit_params(2);
    rf_depths(end+1)      = max_dip;
    valid_filenames{end+1}= rawFiles(i).name;
    snr_array(end+1)      = snr;
    depth_array(end+1)    = max_dip;
    width_array(end+1)    = length(v_idx);
    fine_f = linspace(freq(1), freq(end), 1500).';
    fit_curves_cell{end+1}= [fine_f, env_fun(fine_f) - lor(fit_params, fine_f)];
end

numFiles = length(rf_data_cell);
if numFiles == 0, error('No valid dips found.'); end
max_depth_val = max(rf_depths);
rf_depths = rf_depths(:) / max_depth_val;
rf_centers_Hz = rf_centers_Hz(:);
for i = 1:numFiles, rf_data_cell{i}(:,2) = rf_data_cell{i}(:,2) / max_depth_val; end
mean_rf_Hz = mean(rf_centers_Hz);

%% 4. COMPUTE THE GLOBAL MAPPING
hitran_interp = griddedInterpolant(hitran_freq_Hz, hitran_shape, 'linear', 'none');
D = struct('rf_centers_Hz', rf_centers_Hz, 'rf_depths', rf_depths, ...
           'rf_data_cell', {rf_data_cell}, 'sig_hitran_THz', sig_hitran_THz(:), ...
           'sig_intensity', sig_intensity(:), 'hitran_interp', hitran_interp, ...
           'mean_rf_Hz', mean_rf_Hz);
M = compute_map(P, D);

%% 5. BUILD THE UI
build_ui(P, D, M, struct('rf_raw_cell',{rf_raw_cell}, 'fit_curves_cell',{fit_curves_cell}, ...
    'valid_filenames',{valid_filenames}, 'snr_array',snr_array, 'depth_array',depth_array, ...
    'width_array',width_array, 'hitran_freq_Hz',hitran_freq_Hz, 'hitran_shape',hitran_shape, ...
    'folderPath',folderPath));

% =========================================================================
%% ===========================  HELPER FUNCTIONS  =========================
% =========================================================================

% ----------------------- MAPPING COMPUTATION ----------------------------
function M = compute_map(P, D)
% Band-anchored 2D grid (stretch x centroid) + profile and line (ICP) fits.
    drep = P.f_rep_1 - P.f_rep_2;  dceo = P.f_ceo_1 - P.f_ceo_2;
    m_base = -P.f_rep_1 / drep;
    b_nom  = P.f_ceo_1 - (dceo * P.f_rep_1 / drep);

    rfc = D.rf_centers_Hz(:); rfd = D.rf_depths(:); mean_rf = D.mean_rf_Hz;
    sig = D.sig_hitran_THz(:); sigI = D.sig_intensity(:);
    hit_lo = min(sig); hit_hi = max(sig);

    % line envelope used for the coarse band-lock
    W = 25e9; gf = ((hit_lo-2.5):2e-4:(hit_hi+2.5)).'*1e12; gc = zeros(size(gf));
    sf = sig*1e12;
    for k = 1:numel(sf)
        near = abs(gf - sf(k)) < 5*W;
        gc(near) = gc(near) + sigI(k).*exp(-0.5*((gf(near)-sf(k))/W).^2);
    end
    gc = gc / max(gc);  gi = griddedInterpolant(gf, gc, 'linear', 'none');

    centers = ((hit_lo-0.5):5e-4:(hit_hi+0.5))*1e12;
    rfd_c   = rfc - mean_rf;  rfv = rfd.^3;
    rfc_min = min(rfd_c); rfc_max = max(rfd_c);
    blo = hit_lo*1e12; bhi = hit_hi*1e12; marg = 0.3e12;
    m_shifts = -60000:100:60000;
    best = -inf; m_macro = m_base; m_macro_ppm = 0; center_macro = (hit_lo+hit_hi)/2*1e12;
    for i = 1:numel(m_shifts)
        mt = m_base*(1 + m_shifts(i)*1e-6);
        tf = mt*rfd_c + centers;  ex = gi(tf); ex(isnan(ex)) = 0;
        sc = sum(rfv .* (ex.^3), 1);
        lo_c = blo - marg - mt*rfc_min;  hi_c = bhi + marg - mt*rfc_max;
        sc(centers < lo_c | centers > hi_c) = -inf;
        [ms, mi] = max(sc);
        if ms > best, best = ms; center_macro = centers(mi); m_macro = mt; m_macro_ppm = m_shifts(i); end
    end
    if isinf(best), m_macro = m_base; m_macro_ppm = 0; center_macro = 0.5*(blo+bhi); end
    b_macro = center_macro - m_macro*mean_rf;
    b_base  = b_macro;                                  % re-base to the band

    band_center = mean((m_macro*rfc + b_macro)/1e12);
    lm = (sig > band_center-5) & (sig < band_center+5);
    M.filtered_lines_THz = sig(lm);  M.filtered_intensity = sigI(lm);

    % line fit: resolve the mode by teeth-RMS, then fit curvature of the order
    if ~isfield(P,'map_order'), P.map_order = 1; end
    ord = P.map_order;
    uc = mean_rf;  us = std(rfc);  if ~isfinite(us) || us <= 0, us = 1; end
    [mlf, blf, a2, a3] = line_fit_best(m_base, center_macro, mean_rf, uc, us, ord, rfc, sig, P.f_rep_1, 60000);
    M.m_line_ppm = (mlf/m_base - 1)*1e6;
    M.b_line_THz = (blf - b_base)/1e12;
    M.a2 = a2;  M.a3 = a3;  M.uc = uc;  M.us = us;  M.order = ord;

    M.m_base = m_base;  M.b_base_Hz = b_base;  M.b_base_nominal_Hz = b_nom;
    M.m_macro_ppm = m_macro_ppm;
    M.b_theo_THz  = (b_nom - b_base)/1e12;
end

function [m, b, a2, a3, n_v] = fit_poly(m, b, a2, a3, uc, us, ord, rfc, sig, f_rep, n_iter)
% Iterative closest point with a least-squares POLYNOMIAL of order `ord`
% (1=linear, 2=quadratic, 3=cubic). One-to-one nearest assignment + LS, to
% convergence. The linear part is returned as (m,b); the curvature as (a2,a3),
% the coefficients of u^2 and u^3 with u=(f_rf-uc)/us. This models real comb
% dispersion (the two combs' teeth are not perfectly equidistant), which makes
% f_opt a non-linear function of f_rf.
    rfc = rfc(:); sig = sig(:); frep_THz = f_rep/1e12; u = (rfc-uc)/us; n_v = 0;
    for it = 1:n_iter
        fo = m*rfc + b;
        if ord>=2, fo = fo + a2*u.^2; end
        if ord>=3, fo = fo + a3*u.^3; end
        tgt = assign_lines(fo/1e12, sig, frep_THz);
        g = ~isnan(tgt); n_v = sum(g);
        if n_v < ord+1, return; end
        ug = u(g);  U = ones(n_v, ord+1);
        for kk = 1:ord, U(:,kk+1) = ug.^kk; end       % [1 u u^2 ...]
        coef = U \ (tgt(g)*1e12);                     % a0,a1,a2u,a3u in u-space
        mN = coef(2)/us;  bN = coef(1) - coef(2)*uc/us;
        a2N = 0; a3N = 0;
        if ord>=2, a2N = coef(3); end
        if ord>=3, a3N = coef(4); end
        cv = abs(mN-m) < 1e-7*max(abs(m),1) && abs(bN-b) < 1 && ...
             abs(a2N-a2) < 1 && abs(a3N-a3) < 1;
        m = mN; b = bN; a2 = a2N; a3 = a3N;
        if cv, return; end
    end
end

function [m, b, a2, a3, nv, bestrms] = line_fit_best(m_base, center_macro, mean_rf, uc, us, ord, rfc, sig, f_rep, ppm_range)
% Robust line fit that RESOLVES THE MODE by teeth-RMS (not amplitude overlap):
% scan stretch x integer-line offset, refine each with a RIGID linear ICP
% (order 1 can't overfit), keep the lowest-RMS registration, THEN fit the
% curvature of the requested order at that mode. This reliably finds the
% near-perfect mode that a flat amplitude score is blind to.
    frep_THz = f_rep/1e12; rfc = rfc(:); sig = sig(:);
    dl = median(diff(sort(sig)))*1e12;  if ~isfinite(dl) || dl<=0, dl = 0; end
    ppm_list = -ppm_range:200:ppm_range;
    boff = dl*[-1 -0.5 0 0.5 1];
    need = max(2, round(0.5*numel(rfc)));
    bestrms = inf; m = m_base; b = center_macro - m_base*mean_rf; a2 = 0; a3 = 0; nv = 0;
    for pp = ppm_list
        m0 = m_base*(1 + pp*1e-6);
        for bo = boff
            b0 = center_macro - m0*mean_rf + bo;
            [mm, bb, ~, ~, nn] = fit_poly(m0, b0, 0, 0, uc, us, 1, rfc, sig, f_rep, 40);
            if nn < need, continue; end
            tt = nearest_teeth((mm*rfc + bb)/1e12, sig, frep_THz);
            r = rms(tt);
            if r < bestrms, bestrms = r; m = mm; b = bb; nv = nn; end
        end
    end
    if ord >= 2 && nv >= ord+1            % refine curvature at the chosen mode
        [m, b, a2, a3, nv] = fit_poly(m, b, 0, 0, uc, us, ord, rfc, sig, f_rep, 80);
    end
end

function fo = map_eval(S, frf)
% Apply the current mapping (linear m,b from the sliders + fitted curvature
% a2,a3 of the selected order) to RF frequencies frf [Hz] -> optical [Hz].
    frf = frf(:);
    m = S.m_base*(1+S.sM.Value*1e-6);  b = S.b_base_Hz + S.sB.Value*1e12;
    fo = m*frf + b;
    if S.order>=2
        u = (frf - S.uc)/S.us;
        fo = fo + S.a2*u.^2;
        if S.order>=3, fo = fo + S.a3*u.^3; end
    end
end

function tt = nearest_teeth(mapped_THz, sig, frep_THz)
% Signed distance (in teeth) from each mapped dip to its nearest HITRAN line.
    mp = mapped_THz(:); sig = sig(:); tt = zeros(numel(mp),1);
    for k = 1:numel(mp), [~,j] = min(abs(sig - mp(k))); tt(k) = (mp(k)-sig(j))/frep_THz; end
end

function [bv, envfun] = comb_envelope(tf, tp)
% Robust upper-envelope of the comb power: an asymmetric iterative polynomial
% (deg<=2) through the UNABSORBED teeth. Absorbed teeth sit well below the
% envelope and are rejected each pass, so the fit tracks the comb's Gaussian
% power top (slope + gentle curvature) rather than dipping into the absorption.
    mu = mean(tf); sd = max(std(tf), eps);
    tn = (tf - mu)/sd;  deg = min(2, numel(tp)-1);
    inl = true(size(tp));  pc = polyfit(tn, tp, deg);
    for it = 1:8
        pc = polyfit(tn(inl), tp(inl), deg);
        rr = tp - polyval(pc, tn);
        sg = std(rr(inl));  if sg < eps, break; end
        ni = rr > -1.0*sg;  ni(rr >= 0) = true;        % reject only well-below teeth
        if sum(ni) <= deg+1 || isequal(ni, inl), break; end
        inl = ni;
    end
    bv = polyval(pc, tn);
    envfun = @(f) polyval(pc, (f(:)-mu)/sd);
end

function tgt = assign_lines(f0_THz, sig, frep_THz)
% One-to-one global-nearest-first assignment of dips to distinct HITRAN lines.
    f0 = f0_THz(:); sig = sig(:); nd = numel(f0); nL = numel(sig); HUNT = 0.10;
    cand = zeros(0,3);
    for i = 1:nd
        dd = sig - f0(i); j = find(abs(dd) <= HUNT);
        for jj = j(:).', cand(end+1,:) = [i, jj, abs(dd(jj))/frep_THz]; end
    end
    tgt = nan(nd,1);
    if ~isempty(cand)
        [~,o] = sort(cand(:,3)); du = false(nd,1); lu = false(nL,1);
        for r = o(:).'
            i = cand(r,1); jj = cand(r,2);
            if du(i) || lu(jj), continue; end
            tgt(i) = sig(jj); du(i) = true; lu(jj) = true;
        end
    end
end

function cst = lor_cost(p, f, d, lb, ub)
    if p(1)<lb(1)||p(1)>ub(1)||p(2)<lb(2)||p(2)>ub(2)||abs(p(3))<lb(3)||abs(p(3))>ub(3)
        cst = 1e18; return;
    end
    lorv = p(1).*(p(3)/2).^2 ./ ((f - p(2)).^2 + (p(3)/2).^2);
    cst  = sum((d - lorv).^2);
end

% ----------------------- UI CONSTRUCTION --------------------------------
function build_ui(P, D, M, X)
    nF = numel(X.valid_filenames);
    BG = [0.96 0.96 0.97];

    % ============ MAIN FIGURE ============
    fig = figure('Name','RF -> Optical Mapping', 'NumberTitle','off', ...
        'Color', BG, 'Position', [30 90 1320 780]);

    axC = axes('Parent',fig,'Position',[0.055 0.585 0.62 0.365]); hold(axC,'on'); grid(axC,'on'); box(axC,'on');
    title(axC,'Mapping Comparator','FontWeight','bold');
    ylabel(axC,'Normalised Absorption'); xlabel(axC,'Optical Frequency (THz)');

    axR = axes('Parent',fig,'Position',[0.055 0.365 0.62 0.14]); hold(axR,'on'); grid(axR,'on'); box(axR,'on');
    title(axR,'Per-dip Residual','FontWeight','bold');
    ylabel(axR,'Mismatch (teeth)'); xlabel(axR,'Optical Frequency (THz)');
    yline(axR,0,'k-');

    % HITRAN objects
    h_voigt = plot(axC, X.hitran_freq_Hz/1e12, X.hitran_shape, 'k', 'LineWidth',1.4);
    h_lines = stem(axC, M.filtered_lines_THz, M.filtered_intensity, 'k', 'Marker','.', ...
        'MarkerSize',12, 'LineWidth',1.4, 'Visible','off');

    % RF objects
    pC = gobjects(nF,1); pS = gobjects(nF,1);
    for i = 1:nF
        pC(i) = plot(axC, D.rf_data_cell{i}(:,1), D.rf_data_cell{i}(:,2), 'r', 'LineWidth',2);
        pS(i) = stem(axC, D.rf_centers_Hz(i), D.rf_depths(i), 'Color',[0.85 0 0], ...
            'Marker','o', 'MarkerFaceColor',[0.85 0 0], 'LineWidth',1.5, 'Visible','off');
    end
    % residual scatter
    hRes = scatter(axR, D.rf_centers_Hz, zeros(nF,1), 36, repmat([0.85 0 0],nF,1), 'filled');

    % colour legend on the comparator (traffic-light by teeth mismatch)
    lg1 = plot(axC,nan,nan,'o','MarkerFaceColor',[0.00 0.65 0.10],'MarkerEdgeColor','none','DisplayName','\leq1 tooth');
    lg2 = plot(axC,nan,nan,'o','MarkerFaceColor',[0.92 0.62 0.00],'MarkerEdgeColor','none','DisplayName','\leq5 teeth');
    lg3 = plot(axC,nan,nan,'o','MarkerFaceColor',[0.85 0.05 0.05],'MarkerEdgeColor','none','DisplayName','>5 teeth');
    lgd = legend(axC,[lg1 lg2 lg3],'Location','northwest','FontSize',8,'AutoUpdate','off');
    title(lgd,'mapped-dip mismatch');

    % --------- SLIDERS ----------
    uicontrol(fig,'Style','text','Units','normalized','Position',[0.055 0.255 0.10 0.022], ...
        'String','Stretch  m','FontWeight','bold','HorizontalAlignment','left','BackgroundColor',BG);
    sM = uicontrol(fig,'Style','slider','Units','normalized','Position',[0.155 0.257 0.335 0.024], ...
        'Min',-100000,'Max',100000,'Value',M.m_line_ppm);
    tM = uicontrol(fig,'Style','text','Units','normalized','Position',[0.495 0.255 0.18 0.024], ...
        'String','','HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8);

    uicontrol(fig,'Style','text','Units','normalized','Position',[0.055 0.213 0.10 0.022], ...
        'String','Offset  b','FontWeight','bold','HorizontalAlignment','left','BackgroundColor',BG);
    sB = uicontrol(fig,'Style','slider','Units','normalized','Position',[0.155 0.215 0.335 0.024], ...
        'Min',-30,'Max',30,'Value',M.b_line_THz);
    tB = uicontrol(fig,'Style','text','Units','normalized','Position',[0.495 0.213 0.18 0.024], ...
        'String','','HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8);

    % --------- MAPPING ORDER SELECTOR ----------
    uicontrol(fig,'Style','text','Units','normalized','Position',[0.055 0.168 0.10 0.022], ...
        'String','Map order','FontWeight','bold','HorizontalAlignment','left','BackgroundColor',BG);
    dOrd = uicontrol(fig,'Style','popupmenu','Units','normalized','Position',[0.155 0.169 0.14 0.026], ...
        'String',{'1  (linear)','2  (quadratic)','3  (cubic)'},'Value',M.order, ...
        'TooltipString','Polynomial order of f_opt(f_rf); >1 corrects comb-dispersion curvature');
    uicontrol(fig,'Style','text','Units','normalized','Position',[0.305 0.165 0.37 0.026], ...
        'String','(>1 fits the comb-dispersion curvature; refit on change)', ...
        'HorizontalAlignment','left','BackgroundColor',BG,'FontSize',8,'ForegroundColor',[0.4 0.4 0.4]);

    % --------- ACTION BUTTONS (bottom row) ----------
    bw = 0.099; bx = 0.055; by = 0.045; bh = 0.05; gap = 0.005;
    btnTheo  = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Theoretical', ...
        'Position',[bx by bw bh],'TooltipString','Pure comb-equation mapping (f_{rep}, f_{ceo} only)');
    btnLine  = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Auto Fit', ...
        'Position',[bx+(bw+gap) by bw bh],'TooltipString','Auto-find the best mapping: teeth-RMS mode search + polynomial of the chosen order');
    btnSnap  = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Local Peak Snap', ...
        'Position',[bx+2*(bw+gap) by bw bh],'TooltipString','Re-converge the fit from the current alignment');
    btnMdn   = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Prev Mode', ...
        'Position',[bx+3*(bw+gap) by bw bh],'TooltipString','Shift registration one HITRAN line down, then snap');
    btnMup   = uicontrol(fig,'Style','pushbutton','Units','normalized','String','Next Mode', ...
        'Position',[bx+4*(bw+gap) by bw bh],'TooltipString','Shift registration one HITRAN line up, then snap');

    % --------- RIGHT COLUMN PANELS ----------
    % View  (display style is decoupled from the mapping — both views use the same fit)
    pView = uipanel(fig,'Title','View','FontWeight','bold','BackgroundColor',BG, ...
        'Position',[0.70 0.745 0.285 0.21]);
    uicontrol(pView,'Style','text','Units','normalized','Position',[0.06 0.83 0.28 0.11], ...
        'String','View as:','HorizontalAlignment','left','BackgroundColor',BG);
    dViewAs = uicontrol(pView,'Style','popupmenu','Units','normalized','Position',[0.36 0.84 0.58 0.12], ...
        'String',{'Lines','Profile'},'Value',1, ...
        'TooltipString','Display style only; the mapping is identical either way');
    cbH = uicontrol(pView,'Style','checkbox','Units','normalized','Position',[0.06 0.66 0.9 0.13], ...
        'String','Show HITRAN','Value',1,'BackgroundColor',BG);
    cbR = uicontrol(pView,'Style','checkbox','Units','normalized','Position',[0.06 0.52 0.9 0.13], ...
        'String','Show RF data','Value',1,'BackgroundColor',BG);
    cbCol = uicontrol(pView,'Style','checkbox','Units','normalized','Position',[0.06 0.38 0.9 0.13], ...
        'String','Colour dips by mismatch','Value',1,'BackgroundColor',BG);
    uicontrol(pView,'Style','text','Units','normalized','Position',[0.06 0.20 0.5 0.11], ...
        'String','Highlight file:','HorizontalAlignment','left','BackgroundColor',BG);
    dHl = uicontrol(pView,'Style','popupmenu','Units','normalized','Position',[0.06 0.04 0.88 0.13], ...
        'String',[{'None'}, X.valid_filenames]);

    % Quality
    pQual = uipanel(fig,'Title','Mapping Quality','FontWeight','bold','BackgroundColor',BG, ...
        'Position',[0.70 0.565 0.285 0.17]);
    tQual = uicontrol(pQual,'Style','text','Units','normalized','Position',[0.04 0.04 0.93 0.92], ...
        'String','','HorizontalAlignment','left','BackgroundColor',[1 1 1],'FontName','Menlo','FontSize',9);

    % Mapping Equation (live coefficients of the chosen-order polynomial)
    pEq = uipanel(fig,'Title','Mapping Equation','FontWeight','bold','BackgroundColor',BG, ...
        'Position',[0.70 0.385 0.285 0.17]);
    tEq = uicontrol(pEq,'Style','text','Units','normalized','Position',[0.04 0.04 0.93 0.92], ...
        'String','','HorizontalAlignment','left','BackgroundColor',[1 1 1],'FontName','Menlo','FontSize',9);

    % Comb parameters
    pPar = uipanel(fig,'Title','Comb Parameters (Hz)','FontWeight','bold','BackgroundColor',BG, ...
        'Position',[0.70 0.125 0.285 0.25]);
    lblP = {'f_{rep,1}','f_{ceo,1}','f_{rep,2}','f_{ceo,2}'};
    valP = {P.f_rep_1, P.f_ceo_1, P.f_rep_2, P.f_ceo_2};
    eP = gobjects(4,1);
    for r = 1:4
        yy = 0.80 - (r-1)*0.18;
        uicontrol(pPar,'Style','text','Units','normalized','Position',[0.04 yy 0.30 0.13], ...
            'String',lblP{r},'HorizontalAlignment','left','BackgroundColor',BG);
        eP(r) = uicontrol(pPar,'Style','edit','Units','normalized','Position',[0.36 yy 0.60 0.14], ...
            'String',sprintf('%.6f',valP{r}),'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
    end
    btnRec = uicontrol(pPar,'Style','pushbutton','Units','normalized','Position',[0.10 0.05 0.80 0.16], ...
        'String','Recompute Mapping','FontWeight','bold');

    % Export
    btnExp = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[0.70 0.045 0.285 0.06], ...
        'String','Export Results...','FontWeight','bold');

    % ============ VERIFIER FIGURE ============
    figV = figure('Name','Lorentzian Fit Verifier','NumberTitle','off','Color',BG,'Position',[30 90 760 360]);
    axV = axes('Parent',figV,'Position',[0.08 0.16 0.64 0.74]); hold(axV,'on'); grid(axV,'on'); box(axV,'on');
    ylabel(axV,'Power (dBm)'); xlabel(axV,'RF Freq (MHz)'); title(axV,'Lorentzian Fit Verifier','FontWeight','bold');
    pr = gobjects(nF,1); pt = gobjects(nF,1); pfl = gobjects(nF,1); pff = gobjects(nF,1); pcl = gobjects(nF,1);
    for i = 1:nF
        pr(i) = plot(axV, X.rf_raw_cell{i}(:,1)/1e6, X.rf_raw_cell{i}(:,2), 'Color',[0.6 0.8 1.0], 'Visible','off');
        pt(i) = plot(axV, D.rf_data_cell{i}(:,1)/1e6, D.rf_data_cell{i}(:,3), 'ro', 'MarkerSize',4, 'Visible','off');
        pfl(i)= fill(axV, [D.rf_data_cell{i}(:,1)/1e6; flip(D.rf_data_cell{i}(:,1)/1e6)], ...
                     [D.rf_data_cell{i}(:,4); flip(D.rf_data_cell{i}(:,3))], 'r', 'FaceAlpha',0.2, 'EdgeColor','none','Visible','off');
        pff(i)= plot(axV, X.fit_curves_cell{i}(:,1)/1e6, X.fit_curves_cell{i}(:,2), 'k-', 'LineWidth',2, 'Visible','off');
        pcl(i)= xline(axV, D.rf_centers_Hz(i)/1e6, 'm--', 'LineWidth',1.5, 'Label','fit centre', ...
                     'LabelOrientation','horizontal', 'Visible','off');
    end
    tV = uicontrol(figV,'Style','text','Units','normalized','Position',[0.75 0.30 0.23 0.30], ...
        'String','','FontSize',11,'FontWeight','bold','HorizontalAlignment','left','BackgroundColor',BG);
    dV = uicontrol(figV,'Style','popupmenu','Units','normalized','Position',[0.75 0.74 0.23 0.08], 'String',X.valid_filenames);
    uicontrol(figV,'Style','pushbutton','String','< Prev','Units','normalized','Position',[0.75 0.63 0.10 0.08]);
    uicontrol(figV,'Style','pushbutton','String','Next >','Units','normalized','Position',[0.88 0.63 0.10 0.08]);

    % ============ MISMATCH FIGURE ============
    figM = figure('Name','Individual Dip Mode Mismatch','NumberTitle','off','Color',BG,'Position',[800 90 560 360]);
    t3 = uicontrol(figM,'Style','text','Units','normalized','Position',[0.10 0.92 0.80 0.06], ...
        'FontSize',12,'FontWeight','bold','String','Dip Mismatch','HorizontalAlignment','center','BackgroundColor',BG);
    uicontrol(figM,'Style','text','String','Examine Dip:','Units','normalized','Position',[0.10 0.83 0.18 0.06], ...
        'HorizontalAlignment','left','BackgroundColor',BG);
    d3 = uicontrol(figM,'Style','popupmenu','String',X.valid_filenames,'Units','normalized','Position',[0.27 0.84 0.32 0.05]);
    uicontrol(figM,'Style','pushbutton','String','< Prev','Units','normalized','Position',[0.62 0.83 0.12 0.06]);
    uicontrol(figM,'Style','pushbutton','String','Next >','Units','normalized','Position',[0.76 0.83 0.12 0.06]);
    a3 = axes('Parent',figM,'Position',[0.12 0.14 0.76 0.56]); hold(a3,'on'); grid(a3,'on'); box(a3,'on');
    a3t = axes('Parent',figM,'Position',[0.12 0.14 0.76 0.56],'XAxisLocation','top','YAxisLocation','right', ...
        'Color','none','YTick',[],'HitTest','off','PickableParts','none');
    xlabel(a3,'\Delta n (teeth from target HITRAN)'); ylabel(a3,'Normalised Absorption');
    xlabel(a3t,'Optical Frequency (THz)');

    % --------- ASSEMBLE STATE ----------
    S = struct();
    S.P = P; S.D = D;
    S.m_base = M.m_base; S.b_base_Hz = M.b_base_Hz; S.b_base_nominal_Hz = M.b_base_nominal_Hz;
    S.f_rep_1 = P.f_rep_1; S.mean_rf_Hz = D.mean_rf_Hz;
    S.m_macro_ppm = M.m_macro_ppm;
    S.m_line_ppm = M.m_line_ppm; S.b_line_THz = M.b_line_THz; S.b_theo_THz = M.b_theo_THz;
    S.a2 = M.a2; S.a3 = M.a3; S.uc = M.uc; S.us = M.us; S.order = M.order;
    S.filtered_lines_THz = M.filtered_lines_THz; S.filtered_intensity = M.filtered_intensity;
    S.sig_hitran_THz = D.sig_hitran_THz; S.sig_intensity = D.sig_intensity;
    S.rf_centers_Hz = D.rf_centers_Hz; S.rf_depths = D.rf_depths; S.rf_data_cell = D.rf_data_cell;
    S.valid_filenames = X.valid_filenames; S.folderPath = X.folderPath;
    S.snr_array = X.snr_array; S.depth_array = X.depth_array; S.width_array = X.width_array;
    S.rf_raw_cell = X.rf_raw_cell;
    % handles
    S.fig=fig; S.axC=axC; S.axR=axR; S.h_voigt=h_voigt; S.h_lines=h_lines;
    S.pC=pC; S.pS=pS; S.hRes=hRes; S.sM=sM; S.sB=sB; S.tM=tM; S.tB=tB;
    S.cbH=cbH; S.cbR=cbR; S.cbCol=cbCol; S.dHl=dHl; S.tQual=tQual; S.eP=eP; S.dOrd=dOrd;
    S.dViewAs=dViewAs; S.tEq=tEq;
    S.axV=axV; S.pr=pr; S.pt=pt; S.pfl=pfl; S.pff=pff; S.pcl=pcl; S.tV=tV; S.dV=dV;
    S.figM=figM; S.axM=a3; S.axMt=a3t; S.t3=t3; S.d3=d3;
    S.view='lines'; S.guard=false;
    S.anchor_Hz = (S.m_base*(1+S.m_line_ppm*1e-6))*S.mean_rf_Hz + (S.b_base_Hz + S.b_line_THz*1e12);
    setappdata(fig,'S',S);

    % --------- WIRE CALLBACKS ----------
    btnTheo.Callback = @(s,e) snap_to(fig,'theo');
    btnLine.Callback = @(s,e) snap_to(fig,'auto');
    btnSnap.Callback = @(s,e) local_snap(fig);
    btnMdn.Callback  = @(s,e) mode_step(fig,-1);
    btnMup.Callback  = @(s,e) mode_step(fig,+1);
    btnRec.Callback  = @(s,e) recompute(fig);
    btnExp.Callback  = @(s,e) export_results(fig);
    cbH.Callback   = @(s,e) apply_view(fig);
    cbR.Callback   = @(s,e) apply_view(fig);
    cbCol.Callback = @(s,e) refresh(fig);
    dHl.Callback   = @(s,e) apply_view(fig);
    dViewAs.Callback = @(s,e) apply_view(fig);
    dOrd.Callback  = @(s,e) apply_order(fig);
    addlistener(sM,'Value','PostSet',@(s,e) on_slider(fig,'m'));
    addlistener(sB,'Value','PostSet',@(s,e) on_slider(fig,'b'));
    addlistener(a3,'XLim','PostSet',@(s,e) sync_top_axis(a3,a3t,getappdata(fig,'S')));
    % verifier nav
    set(findobj(figV,'String','< Prev'),'Callback',@(s,e) verif_step(fig,-1));
    set(findobj(figV,'String','Next >'),'Callback',@(s,e) verif_step(fig,+1));
    dV.Callback = @(s,e) verif_step(fig,0);
    % mismatch nav
    set(findobj(figM,'String','< Prev'),'Callback',@(s,e) mism_step(fig,-1));
    set(findobj(figM,'String','Next >'),'Callback',@(s,e) mism_step(fig,+1));
    d3.Callback = @(s,e) refresh(fig);

    % --------- INITIAL DRAW ----------
    apply_view(fig);
    autoscale_comp(fig);
    refresh(fig);
    verif_step(fig,0);
end

% ----------------------- INTERACTION ------------------------------------
function on_slider(fig, which)
    S = getappdata(fig,'S');
    if S.guard, return; end
    if strcmp(which,'m')                                  % stretch pivots about the centroid
        m_new = S.m_base*(1 + S.sM.Value*1e-6);
        b_THz = (S.anchor_Hz - m_new*S.mean_rf_Hz - S.b_base_Hz)/1e12;
        b_THz = max(S.sB.Min, min(S.sB.Max, b_THz));
        S.guard = true; setappdata(fig,'S',S);
        S.sB.Value = b_THz;                               % fires on_slider('b') (guarded)
        S = getappdata(fig,'S'); S.guard = false; setappdata(fig,'S',S);
    else                                                  % offset redefines the centroid anchor
        m_cur = S.m_base*(1 + S.sM.Value*1e-6);
        S.anchor_Hz = m_cur*S.mean_rf_Hz + (S.b_base_Hz + S.sB.Value*1e12);
        setappdata(fig,'S',S);
    end
    refresh(fig);
end

function set_mapping(fig, m_full, b_full, a2, a3, do_autoscale)
    S = getappdata(fig,'S');
    m_ppm = max(S.sM.Min, min(S.sM.Max, (m_full/S.m_base - 1)*1e6));
    b_THz = max(S.sB.Min, min(S.sB.Max, (b_full - S.b_base_Hz)/1e12));
    S.guard = true; S.a2 = a2; S.a3 = a3; setappdata(fig,'S',S);
    S.sM.Value = m_ppm; S.sB.Value = b_THz;
    S = getappdata(fig,'S'); S.guard = false;
    S.anchor_Hz = (S.m_base*(1+m_ppm*1e-6))*S.mean_rf_Hz + (S.b_base_Hz + b_THz*1e12);
    setappdata(fig,'S',S);
    if nargin>5 && do_autoscale, autoscale_comp(fig); end
    refresh(fig);
    drawnow;   % force the slider thumbs to repaint to their new programmatic values
end

function snap_to(fig, kind)
    % Compute a mapping and load it onto the sliders. Display style (lines vs
    % profile) is independent and not touched here.
    S = getappdata(fig,'S');
    switch kind
        case 'theo'                                       % pure comb equations, no curvature
            m = S.m_base; b = S.b_base_nominal_Hz; a2 = 0; a3 = 0;
        case 'auto'                                       % robust mode search + curvature at current order
            m_macro = S.m_base*(1+S.m_macro_ppm*1e-6);
            center_macro = S.b_base_Hz + m_macro*S.mean_rf_Hz;
            [m, b, a2, a3, nv] = line_fit_best(S.m_base, center_macro, S.mean_rf_Hz, S.uc, S.us, ...
                S.order, S.rf_centers_Hz, S.sig_hitran_THz, S.f_rep_1, 60000);
            if nv < S.order+1                             % fall back to stored fit
                m = S.m_base*(1+S.m_line_ppm*1e-6); b = S.b_base_Hz + S.b_line_THz*1e12;
                a2 = S.a2; a3 = S.a3;
            end
    end
    set_mapping(fig, m, b, a2, a3, true);
end

function local_snap(fig)
    S = getappdata(fig,'S');
    m0 = S.m_base*(1+S.sM.Value*1e-6); b0 = S.b_base_Hz + S.sB.Value*1e12;
    frep_THz = S.f_rep_1/1e12;
    tb = nearest_teeth(map_eval(S, S.rf_centers_Hz)/1e12, S.sig_hitran_THz, frep_THz);
    [mN, bN, a2, a3, nv] = fit_poly(m0, b0, S.a2, S.a3, S.uc, S.us, S.order, ...
        S.rf_centers_Hz, S.sig_hitran_THz, S.f_rep_1, 80);
    if nv < S.order+1
        fprintf('[Local Peak Snap] <%d dips within range of a line — align red≈black closer first.\n', S.order+1); return
    end
    set_mapping(fig, mN, bN, a2, a3, false);
    S2 = getappdata(fig,'S');
    ta = nearest_teeth(map_eval(S2, S2.rf_centers_Hz)/1e12, S2.sig_hitran_THz, frep_THz);
    fprintf('[Local Peak Snap] order %d, converged on %d dips: RMS %.2f -> %.2f teeth\n', S.order, nv, rms(tb), rms(ta));
end

function mode_step(fig, dir)
    % Shift the registration by ~one HITRAN line spacing, then snap (poly ICP).
    S = getappdata(fig,'S');
    m0 = S.m_base*(1+S.sM.Value*1e-6); b0 = S.b_base_Hz + S.sB.Value*1e12;
    dl = median(diff(sort(S.filtered_lines_THz)))*1e12;     % Hz, one line spacing
    if ~isfinite(dl) || dl<=0, dl = median(diff(sort(S.sig_hitran_THz)))*1e12; end
    [mN, bN, a2, a3, nv] = fit_poly(m0, b0 + dir*dl, S.a2, S.a3, S.uc, S.us, S.order, ...
        S.rf_centers_Hz, S.sig_hitran_THz, S.f_rep_1, 80);
    if nv < S.order+1, fprintf('[Mode step] no lines in range after shift.\n'); return; end
    set_mapping(fig, mN, bN, a2, a3, false);
    fprintf('[Mode step] shifted %+d line; snapped (order %d).\n', dir, S.order);
end

function apply_order(fig)
    % Re-fit the mapping polynomial at the newly selected order, keeping the
    % current mode (linear position); the curvature terms are recomputed.
    S = getappdata(fig,'S');
    S.order = S.dOrd.Value; setappdata(fig,'S',S);       % 1,2,3
    m0 = S.m_base*(1+S.sM.Value*1e-6); b0 = S.b_base_Hz + S.sB.Value*1e12;
    [m, b, a2, a3, nv] = fit_poly(m0, b0, 0, 0, S.uc, S.us, S.order, ...
        S.rf_centers_Hz, S.sig_hitran_THz, S.f_rep_1, 80);
    if nv >= S.order+1, set_mapping(fig, m, b, a2, a3, false);
    else,               set_mapping(fig, m0, b0, 0, 0, false); end
    fprintf('[Map order] set to %d.\n', S.order);
end

function recompute(fig)
    S = getappdata(fig,'S');
    P = S.P;
    P.f_rep_1 = str2double(get(S.eP(1),'String'));
    P.f_ceo_1 = str2double(get(S.eP(2),'String'));
    P.f_rep_2 = str2double(get(S.eP(3),'String'));
    P.f_ceo_2 = str2double(get(S.eP(4),'String'));
    if any(~isfinite([P.f_rep_1 P.f_ceo_1 P.f_rep_2 P.f_ceo_2])) || P.f_rep_1==P.f_rep_2
        warndlg('Invalid comb parameters (check f_{rep,1} \neq f_{rep,2}).','Recompute'); return;
    end
    fprintf('Recomputing mapping with edited comb parameters...\n');
    P.map_order = S.order;                                % keep the user's chosen order
    M = compute_map(P, S.D);
    S.P = P; S.f_rep_1 = P.f_rep_1;
    S.m_base=M.m_base; S.b_base_Hz=M.b_base_Hz; S.b_base_nominal_Hz=M.b_base_nominal_Hz;
    S.m_macro_ppm=M.m_macro_ppm;
    S.m_line_ppm=M.m_line_ppm; S.b_line_THz=M.b_line_THz; S.b_theo_THz=M.b_theo_THz;
    S.a2=M.a2; S.a3=M.a3; S.uc=M.uc; S.us=M.us; S.order=M.order;
    S.filtered_lines_THz=M.filtered_lines_THz; S.filtered_intensity=M.filtered_intensity;
    setappdata(fig,'S',S);
    set(S.h_lines,'XData',M.filtered_lines_THz,'YData',M.filtered_intensity);
    snap_to(fig,'auto');     % jump to the new best fit
end

% ----------------------- VIEW & REFRESH ---------------------------------
function apply_view(fig)
    % Set visibility of HITRAN + RF objects from the "View as" popup, the
    % show/hide checkboxes, and the highlight selection. Display style ONLY —
    % the mapping is unchanged.
    S = getappdata(fig,'S');
    S.view = tern(S.dViewAs.Value==1,'lines','profile'); setappdata(fig,'S',S);
    isLine = strcmp(S.view,'lines');
    if S.cbH.Value
        set(S.h_voigt,'Visible', tern(isLine,'off','on'));
        set(S.h_lines,'Visible', tern(isLine,'on','off'));
    else
        set(S.h_voigt,'Visible','off'); set(S.h_lines,'Visible','off');
    end
    set(S.pC,'Visible','off'); set(S.pS,'Visible','off');
    if S.cbR.Value
        hl = S.dHl.Value;                                   % 1 = None
        if hl <= 1
            if isLine, set(S.pS,'Visible','on'); else, set(S.pC,'Visible','on'); end
        else
            if isLine, set(S.pS(hl-1),'Visible','on'); else, set(S.pC(hl-1),'Visible','on'); end
        end
    end
    refresh(fig);
end

function refresh(fig)
    S = getappdata(fig,'S');
    m = S.m_base*(1+S.sM.Value*1e-6); b = S.b_base_Hz + S.sB.Value*1e12;
    frep_THz = S.f_rep_1/1e12; rfc = S.rf_centers_Hz; rfd = S.rf_depths;
    % label shows the linear part (m used in f_opt = m*f_rf + b + curvature),
    % plus the ppm deviation from m_base (the quantity the slider tunes).
    set(S.tM,'String',sprintf('m = %.3f  (%+.0f ppm)', m, S.sM.Value));
    set(S.tB,'String',sprintf('b = %.4e Hz', b));

    mapped = map_eval(S, rfc)/1e12;
    tt = nearest_teeth(mapped, S.sig_hitran_THz, frep_THz);
    cols = zeros(numel(rfc),3);
    for i = 1:numel(rfc), cols(i,:) = mismatch_color(tt(i)); end
    colOn  = logical(S.cbCol.Value);
    useCol = colOn;                                  % colour applies to both views now
    hl = S.dHl.Value;

    for i = 1:numel(rfc)
        set(S.pC(i),'XData', map_eval(S, S.rf_data_cell{i}(:,1))/1e12);
        set(S.pS(i),'XData',mapped(i),'YData',rfd(i));
        if useCol, c = cols(i,:); else, c = [0.85 0 0]; end
        set(S.pS(i),'Color',c,'MarkerFaceColor',c);
        set(S.pC(i),'Color',c);
    end
    if hl > 1   % highlighted file: blue + bold
        set(S.pS(hl-1),'Color',[0 0.2 0.9],'MarkerFaceColor',[0 0.2 0.9],'LineWidth',2.2);
        set(S.pC(hl-1),'Color',[0 0.2 0.9],'LineWidth',3);
    end

    % residual scatter — also obeys the "colour by mismatch" checkbox
    if colOn, residCols = cols; else, residCols = repmat([0.20 0.40 0.90],numel(rfc),1); end
    set(S.hRes,'XData',mapped,'YData',tt,'CData',residCols);
    mx = max([2; abs(tt(:))]);
    ylim(S.axR,[-1.15*mx 1.15*mx]); xlim(S.axR, xlim(S.axC));

    % curvature contribution across the dips, in teeth (0 for order 1)
    u = (rfc - S.uc)/S.us; hi = zeros(size(u));
    if S.order>=2, hi = hi + S.a2*u.^2; end
    if S.order>=3, hi = hi + S.a3*u.^3; end
    curv_teeth = (max(hi)-min(hi))/S.f_rep_1;
    ordname = {'linear','quadratic','cubic'};

    % quality readout
    nmatch = sum(~isnan(assign_lines(mapped, S.sig_hitran_THz, frep_THz)));
    S.tQual.String = sprintf(['RMS mismatch  : %6.2f teeth\nMean |mismatch|: %6.2f teeth\n' ...
        'Max mismatch  : %6.2f teeth\nMatched 1:1   : %4d / %d\n\n' ...
        'model     : order %d (%s)\ncurv span : %7.1f teeth\nstretch m : %+.0f ppm\noffset  b : %.4e Hz'], ...
        rms(tt), mean(abs(tt)), max(abs(tt)), nmatch, numel(rfc), ...
        S.order, ordname{S.order}, curv_teeth, S.sM.Value, b);

    % live mapping equation (absolute-Hz coefficients of the chosen order)
    set(S.tEq,'String', equation_str(S, m, b));

    draw_mismatch(fig);
end

function [c0,c1,c2,c3] = poly_coeffs(S, m, b)
% Expand the mapping (slider linear part m,b + curvature a2,a3 in u-space,
% u=(f_rf-uc)/us) to absolute-frequency coefficients:
%   f_opt = c3*f_rf^3 + c2*f_rf^2 + c1*f_rf + c0   (Hz)
    a2 = S.a2; a3 = S.a3; uc = S.uc; us = S.us; ord = S.order;
    if ord < 3, a3 = 0; end
    if ord < 2, a2 = 0; end
    c3 = a3/us^3;
    c2 = a2/us^2 - 3*a3*uc/us^3;
    c1 = m - 2*a2*uc/us^2 + 3*a3*uc^2/us^3;
    c0 = b + a2*uc^2/us^2 - a3*uc^3/us^3;
end

function str = equation_str(S, m, b)
% Multi-line display of the chosen-order mapping polynomial + its coefficients.
    [c0,c1,c2,c3] = poly_coeffs(S, m, b);
    forms = {'c1*f_rf + c0', ...
             'c2*f_rf^2 + c1*f_rf + c0', ...
             'c3*f_rf^3 + c2*f_rf^2 + c1*f_rf + c0'};
    str = sprintf('f_opt = %s\n(f_rf, f_opt in Hz)\n\n', forms{S.order});
    if S.order>=3, str = [str sprintf('c3 = %+ .8e\n', c3)]; end
    if S.order>=2, str = [str sprintf('c2 = %+ .8e\n', c2)]; end
    str = [str sprintf('c1 = %+ .8e\n', c1)];
    str = [str sprintf('c0 = %+ .8e', c0)];
end

function draw_mismatch(fig)
    S = getappdata(fig,'S'); idx = S.d3.Value; frep = S.f_rep_1;
    rdata = S.rf_data_cell{idx};
    opt_center = map_eval(S, S.rf_centers_Hz(idx))/1e12;
    [~,j] = min(abs(S.sig_hitran_THz - opt_center)); f_h = S.sig_hitran_THz(j);
    S.axM.UserData = f_h;
    dn = (opt_center - f_h)*1e12/frep;
    rel = (map_eval(S, rdata(:,1))/1e12 - f_h)*1e12/frep;
    keepXY = ~isempty(S.axM.Children) && isfield(S,'mism_idx') && isequal(S.mism_idx, idx);
    if keepXY, xl = S.axM.XLim; yl = S.axM.YLim;
    else
        sp = [min([rel(:);0]) max([rel(:);0])]; pad = max(6,0.3*diff(sp));
        xl = [sp(1)-pad sp(2)+pad]; yl = [0 max(rdata(:,2))*1.25+1e-3];
    end
    cla(S.axM); hold(S.axM,'on');
    stem(S.axM, rel, rdata(:,2),'b','Marker','.','LineWidth',1.5,'DisplayName','Mapped comb teeth');
    xline(S.axM,(opt_center-f_h)*1e12/frep,'b--','LineWidth',1.5,'DisplayName',sprintf('Mapped centre (%.6f THz)',opt_center));
    stem(S.axM,0,max(rdata(:,2))*1.1,'r','Marker','diamond','MarkerFaceColor','r','LineWidth',2, ...
        'DisplayName',sprintf('HITRAN target (%.6f THz)',f_h));
    S.t3.String = sprintf('%s  |  Mismatch: %+.2f teeth', S.valid_filenames{idx}, dn);
    legend(S.axM,'Location','northeast');
    xlim(S.axM,xl); ylim(S.axM,yl);
    S.mism_idx = idx; setappdata(fig,'S',S);
    sync_top_axis(S.axM, S.axMt, S);
end

function mism_step(fig, dir)
    S = getappdata(fig,'S'); nv = S.d3.Value+dir;
    if nv>=1 && nv<=numel(S.d3.String), S.d3.Value=nv; setappdata(fig,'S',S); refresh(fig); end
end

function sync_top_axis(a3, a3t, S)
    if ~isempty(a3.UserData)
        a3t.XLim = a3.UserData + (a3.XLim * S.f_rep_1 / 1e12);
    end
end

function autoscale_comp(fig)
    S = getappdata(fig,'S');
    allx = [];
    for i=1:numel(S.rf_centers_Hz), allx = [allx; map_eval(S, S.rf_data_cell{i}(:,1))/1e12]; end
    if ~isempty(allx), xlim(S.axC,[min(allx)-0.1 max(allx)+0.1]); end
    ylim(S.axC,[-0.05 1.15]);
end

% ----------------------- VERIFIER ---------------------------------------
function verif_step(fig, dir)
    S = getappdata(fig,'S'); d = S.dV; nv = d.Value+dir;
    if nv>=1 && nv<=numel(d.String), d.Value=nv; end
    i = d.Value;
    set(S.pr,'Visible','off'); set(S.pt,'Visible','off'); set(S.pfl,'Visible','off');
    set(S.pff,'Visible','off'); set(S.pcl,'Visible','off');
    set([S.pr(i) S.pt(i) S.pfl(i) S.pff(i) S.pcl(i)],'Visible','on');
    xlim(S.axV,[min(S.rf_raw_cell{i}(:,1))/1e6 max(S.rf_raw_cell{i}(:,1))/1e6]);
    S.tV.String = sprintf('SNR: %.1f\n\nDepth: %.1f dB\n\nWidth: %d teeth\n\nCentre:\n%.6f MHz', ...
        S.snr_array(i), S.depth_array(i), S.width_array(i), S.rf_centers_Hz(i)/1e6);
end

% ----------------------- EXPORT -----------------------------------------
function export_results(fig)
    S = getappdata(fig,'S');
    m = S.m_base*(1+S.sM.Value*1e-6); b = S.b_base_Hz + S.sB.Value*1e12;
    frep_THz = S.f_rep_1/1e12;
    mapped = map_eval(S, S.rf_centers_Hz)/1e12;
    tgt = assign_lines(mapped, S.sig_hitran_THz, frep_THz);
    tt  = nearest_teeth(mapped, S.sig_hitran_THz, frep_THz);

    % Expand the mapping to absolute-frequency polynomial coefficients so it can
    % be applied directly:  f_opt = c3*f_rf^3 + c2*f_rf^2 + c1*f_rf + c0  (Hz).
    [c0,c1,c2,c3] = poly_coeffs(S, m, b);
    ord = S.order;
    ordname = {'linear','quadratic','cubic'};

    [fn,pn] = uiputfile('*.csv','Export mapping results', fullfile(S.folderPath,'rf_to_optical_results.csv'));
    if isequal(fn,0), return; end
    fid = fopen(fullfile(pn,fn),'wt');
    fprintf(fid,'# RF -> Optical mapping results\n');
    fprintf(fid,'# Mapping model:  f_opt = c3*f_rf^3 + c2*f_rf^2 + c1*f_rf + c0   (Hz)\n');
    fprintf(fid,'#   polynomial order = %d (%s)\n', ord, ordname{ord});
    fprintf(fid,'#\n');
    fprintf(fid,'# Applied mapping coefficients (f_rf in Hz, f_opt in Hz)\n');
    fprintf(fid,'#   c0 = %.12g\n', c0);
    fprintf(fid,'#   c1 = %.12g\n', c1);
    fprintf(fid,'#   c2 = %.12g   (0 if order<2)\n', c2);
    fprintf(fid,'#   c3 = %.12g   (0 if order<3)\n', c3);
    fprintf(fid,'#\n');
    fprintf(fid,'# Linear part (sliders): f_opt ~ m*f_rf + b (+ curvature for order>1)\n');
    fprintf(fid,'#   m (stretch) = %.10g\n', m);
    fprintf(fid,'#   b (offset)  = %.10g Hz\n', b);
    fprintf(fid,'#\n');
    fprintf(fid,'# Theoretical linear mapping (from comb parameters, no dispersion)\n');
    fprintf(fid,'#   m_theoretical = %.10g\n', S.m_base);
    fprintf(fid,'#   b_theoretical = %.10g Hz\n', S.b_base_nominal_Hz);
    fprintf(fid,'#\n');
    fprintf(fid,'# Comb parameters (Hz)\n');
    fprintf(fid,'#   f_rep_1 = %.6f\n', S.P.f_rep_1);
    fprintf(fid,'#   f_ceo_1 = %.6f\n', S.P.f_ceo_1);
    fprintf(fid,'#   f_rep_2 = %.6f\n', S.P.f_rep_2);
    fprintf(fid,'#   f_ceo_2 = %.6f\n', S.P.f_ceo_2);
    fprintf(fid,'#\n');
    fprintf(fid,'# Mapping quality\n');
    fprintf(fid,'#   RMS mismatch    = %.3f teeth\n', rms(tt));
    fprintf(fid,'#   Mean |mismatch| = %.3f teeth\n', mean(abs(tt)));
    fprintf(fid,'#   Max mismatch    = %.3f teeth\n', max(abs(tt)));
    fprintf(fid,'#   Matched 1:1     = %d / %d\n', sum(~isnan(tgt)), numel(tt));
    fprintf(fid,'#\n');
    fprintf(fid,'file,rf_center_Hz,rf_center_MHz,mapped_THz,assigned_line_THz,mismatch_teeth,depth_dB\n');
    for i = 1:numel(S.rf_centers_Hz)
        if isnan(tgt(i)), tl = NaN; else, tl = tgt(i); end
        fprintf(fid,'%s,%.6f,%.6f,%.9f,%.9f,%.4f,%.3f\n', S.valid_filenames{i}, ...
            S.rf_centers_Hz(i), S.rf_centers_Hz(i)/1e6, mapped(i), tl, tt(i), S.depth_array(i));
    end
    fclose(fid);
    fprintf('Exported %d dips to %s\n', numel(S.rf_centers_Hz), fullfile(pn,fn));
    msgbox(sprintf('Exported %d dips to\n%s', numel(S.rf_centers_Hz), fullfile(pn,fn)),'Export complete');
end

% ----------------------- SMALL UTILITIES --------------------------------
function c = mismatch_color(t)
    a = abs(t);
    if     a < 1, c = [0.00 0.65 0.10];     % green: on the line
    elseif a < 5, c = [0.92 0.62 0.00];     % amber: within a few teeth
    else,         c = [0.85 0.05 0.05];     % red: off
    end
end

function v = tern(cond,a,b)
    if cond, v=a; else, v=b; end
end
