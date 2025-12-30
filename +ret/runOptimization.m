function [S, success] = runOptimization(S)
%   OPTIMIZATION STRATEGY:
%   1. Iterate Pin Diameter (Discrete)
%   2. Iterate Wall Thickness (Outer loop) -> Check Hoop/Axial Stress immediately (Pruning)
%   3. Iterate N_Rows (Integer)
%   4. Iterate N_Pins (Integer) -> Check Bearing/PinShear/Packing (Pruning)
%   5. Iterate Row Spacing (Fine Grid) -> Check Net Tension (Pruning)
%   6. Solve for 'First Row Offset' (Z) analytically to find the exact minimum Z 
%      required to pass Shear-Out. This eliminates the Z loop and guarantees 
%      lowest mass for that combination.

    %% 1. CONFIGURATION & SETUP
    success = false;
    
    % --- SEARCH RESOLUTION (Step Sizes) ---
    % Smaller steps = more precision but slower.
    STEPS.t          = 0.0001;   % Wall thickness step (in)
    STEPS.rowSpacing = 0.0001;   % Row spacing step (in)
    
    % --- EXTRACT CONSTANTS ---
    % Pulling values out of the struct 'S' speeds up loop access
    ID          = S.ID;
    MEOP        = S.MEOP_psi;
    DF_casing   = S.DF_casing;
    DF_pin      = S.DF_pin;
    L_casing    = S.L_casing;
    
    dens_CF     = S.density_CF;
    dens_Al     = S.density_Al;
    dens_Pin    = S.density_pin;
    retRingThk  = S.retRingThk;
    
    % Design Factors & Loads
    p_des_casing = MEOP * DF_casing;
    p_des_pin    = MEOP * DF_pin;
    
    A_bore       = pi * (ID/2)^2;
    F_ax_casing  = p_des_casing * A_bore; % lbf
    F_ax_pin     = p_des_pin    * A_bore; % lbf
    
    % Geometric Factors
    minCircPitchFactor  = S.minCircPitchFactor;
    minAxPitchFactor    = S.minAxialPitchFactor;
    
    % Constraints
    Tg = S.targets; % shearOut_max, netTension_max, etc.
    Bounds = S.optBounds;
    
    % --- MARGIN & PHYSICS CORRECTIONS ---
    % 1. Respect the "Soft Pass" margin from the GUI
    if isfield(S, 'marginFrac')
        mf = S.marginFrac;
    else
        mf = 0.0; 
    end
    
    % 2. Match Pin Strength to GUI (updateAll.m uses 75 KSI, legacy used 43.5)
    PinStrength = 75.0; % KSI
    
    %% 2. INITIALIZE GRID RANGES
    
    % Pin Diameters (Discrete list)
    allowedPins = S.allowedPinDias;
    % Filter pins by bounds
    allowedPins = allowedPins(allowedPins >= Bounds.pinDia(1) & allowedPins <= Bounds.pinDia(2));
    if isempty(allowedPins), error('No pin diameters allowed by current bounds.'); end
    
    % Thickness Range
    t_vals = Bounds.t(1) : STEPS.t : Bounds.t(2);
    
    % Row Counts
    rows_vals = Bounds.nRows(1) : 1 : Bounds.nRows(2);
    
    % Pins Per Row
    ppr_vals  = Bounds.nPinsPerRow(1) : 1 : Bounds.nPinsPerRow(2);
    
    % Row Spacing Range
    sp_vals   = Bounds.rowSpacing(1) : STEPS.rowSpacing : Bounds.rowSpacing(2);
    
    % Best Result Holder
    best.obj  = inf;
    best.S    = S;
    best.mass = inf;
    
    % Progress Bar
    total_iters_est = length(allowedPins) * length(t_vals) * length(rows_vals) * length(ppr_vals) * 0.1; % rough heuristic
    hWait = waitbar(0, 'Running Smart Global Optimization...');
    ctr = 0;
    
    %% 3. EXECUTE SEARCH
    
    % Loop 1: Pin Diameter (Discrete)
    for pIdx = 1:length(allowedPins)
        pinDia = allowedPins(pIdx);
        r_pin  = pinDia/2;
        A_pin  = pi * r_pin^2;
        
        % Pre-calc Pitch Limit
        minAxSpacing = minAxPitchFactor * pinDia;
        
        % Loop 2: Wall Thickness (t)
        for t = t_vals
            
            % --- PRUNING CHECK A: THICKNESS DEPENDENT STRESSES ---
            % Hoop and Axial stresses depend ONLY on t (and ID/Pressure).
            % If t is too thin, NO combination of pins will save it.
            
            OD = ID + 2*t;
            
            % Hoop Stress (Lame Equation)
            % Sigma_h = P * (ro^2 + ri^2)/(ro^2 - ri^2)
            term = ((ID+2*t)^2 + ID^2) / ((ID+2*t)^2 - ID^2);
            Hoop = (p_des_casing * term) / 1000; % KSI
            Axial = Hoop / 2;
            
            % Apply Margin (Allow slight violation if GUI allows it)
            if Hoop > Tg.hoop_max * (1 + mf) || Axial > Tg.axial_max * (1 + mf)
                continue; % Skip all inner loops for this thickness
            end
            
            % Derived Geometry for Packing
            circumference = pi * (ID + t); 
            % Replicating original code: r_o = r_i + t -> circumference = 2*pi*r_o
            r_o = (ID/2) + t;
            circ_outer = 2*pi*r_o; 
            minCircPitch = minCircPitchFactor * pinDia;
            maxPinsGeom  = floor(circ_outer / minCircPitch);
            
            % Update Waitbar occasionally
            ctr = ctr + 1;
            if mod(ctr, 50) == 0
                waitbar(ctr/total_iters_est, hWait, sprintf('Opt: Dia %.3f | t %.3f | Best: %.2f lb', pinDia, t, best.mass));
            end

            % Loop 3: Number of Rows
            for nRows = rows_vals
                
                % Loop 4: Pins Per Row
                for nPPR = ppr_vals
                    
                    % --- PRUNING CHECK B: PACKING ---
                    if nPPR > maxPinsGeom
                        break; % Assuming ppr_vals is sorted ascending, larger values will also fail
                    end
                    
                    pins_oneEnd = nRows * nPPR;
                    if pins_oneEnd == 0, continue; end
                    
                    % --- PRUNING CHECK C: PIN/BEARING STRESSES ---
                    % These depend on count and t, but NOT spacing or Z.
                    
                    F_perPin_casing = F_ax_casing / pins_oneEnd;
                    F_perPin_pin    = F_ax_pin    / pins_oneEnd;
                    
                    % Bearing Stress
                    Bearing = (F_perPin_casing / (pinDia * t)) / 1000;
                    if Bearing > Tg.bearing_max * (1 + mf), continue; end
                    
                    % Pin Shear Stress
                    PinShear = (F_perPin_pin / A_pin) / 1000;
                    if PinShear > Tg.pinShear_max * (1 + mf), continue; end
                    
                    % Pin FOS
                    if PinShear <= 0
                        FOS = inf;
                    else
                        FOS = (PinStrength / PinShear) - 2;
                    end
                    % Min FOS check (allow slight under-performance)
                    if FOS < Tg.pinShearFOS_min * (1 - mf), continue; end
                    
                    % --- LOOP 5: ROW SPACING ---
                    for rsp = sp_vals
                        
                        % Geo Check: Axial Pitch
                        if rsp < minAxSpacing, continue; end
                        
                        % Interaction Factor Logic for Net Tension
                        % NOTE: Must match updateAll.m logic.
                        kInteract = 0.0; 
                        
                        if rsp <= kInteract * pinDia
                            N_eff = nRows;
                        else
                            N_eff = 1;
                        end
                        
                        % Net Tension Stress
                        A_wall_gross    = (pi/4) * (OD^2 - ID^2);
                        A_holes_per_row = pinDia * nPPR * t;
                        A_tens_net      = A_wall_gross - N_eff * A_holes_per_row;
                        
                        if A_tens_net <= 0, A_tens_net = 1e-6; end
                        NetTension = (F_ax_casing / A_tens_net) / 1000;
                        
                        if NetTension > Tg.netTension_max * (1 + mf), continue; end
                        
                        % --- ANALYTIC SOLUTION FOR FIRST ROW OFFSET (Z) ---
                        % Instead of looping Z, we find the MINIMUM Z required to pass Shear Out.
                        % A_shear = t * 2 * nPPR * SUM(spacing*(n-1) + Z) for n=1:Rows
                        
                        % Adjusted Target for Margin
                        target_ShearOut = Tg.shearOut_max * (1 + mf);
                        req_A_shear = (F_ax_casing / 1000) / target_ShearOut;
                        
                        % Term B: Contribution from spacing
                        term_spacing = rsp * (nRows * (nRows-1) / 2);
                        
                        % Solve: t * 2 * nPPR * (Rows*Z + term_spacing) >= req_A_shear
                        % Rows*Z + term_spacing >= req_A_shear / (2*t*nPPR)
                        % Z >= [ (req_A_shear / (2*t*nPPR)) - term_spacing ] / Rows
                        
                        val_RHS = req_A_shear / (2 * t * nPPR);
                        min_Z_shear = (val_RHS - term_spacing) / nRows;
                        
                        % The actual Z must be at least this, AND within bounds
                        final_Z = max(Bounds.firstRowZ(1), min_Z_shear);
                        
                        % Check Upper Bound of Z
                        if final_Z > Bounds.firstRowZ(2)
                            continue; % Impossible to satisfy ShearOut within Z max bound
                        end
                        
                        % Check Length Constraint
                        % The GUI enforces symmetric margins: (firstRowZ + Length + firstRowZ) <= L_casing
                        % Total Occupied = 2*Z + (Rows-1)*Spacing
                        total_occupied_len = 2*final_Z + (nRows-1)*rsp;
                        
                        if total_occupied_len > L_casing
                            continue; % Hardware too long
                        end
                        
                        % --- IF WE REACHED HERE, IT IS A VALID DESIGN ---
                        % Calculate Mass
                        
                        % 1. Casing Mass
                        V_casing_gross = A_wall_gross * S.casingLength;
                        V_holes_casing = (pins_oneEnd*2) * pi * r_pin^2 * t; % *2 for both ends
                        M_casing = (V_casing_gross - V_holes_casing) * dens_CF;
                        
                        % 2. Retention Rings Mass (Al)
                        retLen = total_occupied_len; % Symmetric Z
                        ringOD = ID;
                        ringID = ID - 2*retRingThk;
                        A_ring = (pi/4)*(ringOD^2 - ringID^2);
                        V_ring = A_ring * retLen;
                        V_ring_holes = pins_oneEnd * pi * r_pin^2 * retRingThk;
                        M_rings = (V_ring - V_ring_holes) * dens_Al * 2; % *2 ends
                        
                        % 3. Pins Mass
                        pinLen = t + retRingThk;
                        V_pins = (pins_oneEnd*2) * pi * r_pin^2 * pinLen;
                        M_pins = V_pins * dens_Pin;
                        
                        Mass_total = M_casing + M_rings + M_pins;
                        
                        % Objective Function (Mass + tiny penalty for length to break ties)
                        obj = 1.0 * Mass_total + 0.05 * total_occupied_len;
                        
                        if obj < best.obj
                            best.obj = obj;
                            best.mass = Mass_total;
                            
                            % Store solution parameters
                            best.S.t = t;
                            best.S.nRows = nRows;
                            best.S.nPinsPerRow = nPPR;
                            best.S.rowSpacing = rsp;
                            best.S.firstRowZ = final_Z;
                            best.S.pinDia = pinDia;
                            best.S.pinLen = pinLen;
                            
                            success = true;
                        end
                        
                    end % Spacing
                end % Pins
            end % Rows
        end % Thickness
    end % PinDia
    
    close(hWait);
    
    %% 4. APPLY RESULTS
    if success
        S = best.S;
        
        fprintf('\n=== GLOBAL OPTIMIZATION SUCCESS ===\n');
        fprintf('Total Mass:      %.4f lb\n', best.mass);
        fprintf('Pin Diameter:    %.4f in\n', S.pinDia);
        fprintf('Wall Thickness:  %.4f in\n', S.t);
        fprintf('Configuration:   %d Rows x %d Pins\n', S.nRows, S.nPinsPerRow);
        fprintf('Row Spacing:     %.4f in\n', S.rowSpacing);
        fprintf('First Row Offset:%.4f in\n', S.firstRowZ);
        fprintf('-----------------------------------\n');
        
    else
        fprintf('\n=== OPTIMIZATION FAILED ===\n');
        fprintf('No valid solution found within the specified bounds.\n');
        fprintf('Try widening the bounds (e.g., Min Thickness, Max Rows).\n');
    end

end