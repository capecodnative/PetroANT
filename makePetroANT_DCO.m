% Back-fill elemental *_ppm columns from oxide wt% columns in table T.
% Existing ppm values are kept if they are nonzero and non-NaN.
% Missing/zero ppm values are filled from oxide wt% columns.
%
% Assumes:
%   - oxide columns are lowercase (e.g. sio2, tio2, fe2o3_tot)
%   - ppm columns are lowercase (e.g. si_ppm, ti_ppm, fe_ppm, fetot_ppm)
%   - MolMass() is available, e.g. MolMass('Si'), MolMass('TiO2')

%% Oxide -> element map
map = {
    'sio2',    'si_ppm',  'Si', 1, 'SiO2'
    'tio2',    'ti_ppm',  'Ti', 1, 'TiO2'
    'al2o3',   'al_ppm',  'Al', 2, 'Al2O3'
    'cr2o3',   'cr_ppm',  'Cr', 2, 'Cr2O3'
    'mgo',     'mg_ppm',  'Mg', 1, 'MgO'
    'cao',     'ca_ppm',  'Ca', 1, 'CaO'
    'mno',     'mn_ppm',  'Mn', 1, 'MnO'
    'nio',     'ni_ppm',  'Ni', 1, 'NiO'
    'k2o',     'k_ppm',   'K',  2, 'K2O'
    'na2o',    'na_ppm',  'Na', 2, 'Na2O'
    'sro',     'sr_ppm',  'Sr', 1, 'SrO'
    'p2o5',    'p_ppm',   'P',  2, 'P2O5'
    'co2',     'c_ppm',   'C',  1, 'CO2'
    'so3',     's_ppm',   'S',  1, 'SO3'
    'bao',     'ba_ppm',  'Ba', 1, 'BaO'
};

PetroANT_DCO = PetroANT;

%% Fill standard elemental ppm columns
for i = 1:size(map,1)
    oxideVar   = map{i,1};
    ppmVar     = map{i,2};
    elemSym    = map{i,3};
    nElem      = map{i,4};
    oxideFormula = map{i,5};

    if ~ismember(oxideVar, PetroANT_DCO.Properties.VariableNames)
        continue
    end

    if ~ismember(ppmVar, PetroANT_DCO.Properties.VariableNames)
        PetroANT_DCO.(ppmVar) = nan(height(PetroANT_DCO),1);
    end

    oxide = PetroANT_DCO.(oxideVar);
    ppm   = PetroANT_DCO.(ppmVar);

    ppm_from_oxide = oxide .* 1e4 .* (nElem * MolMass(elemSym) / MolMass(oxideFormula));

    fillIdx = (isnan(ppm) | ppm == 0) & ~isnan(oxide) & oxide ~= 0;
    ppm(fillIdx) = ppm_from_oxide(fillIdx);

    PetroANT_DCO.(ppmVar) = ppm;
end

%% fe_ppm = fe from fe2o3 + feo
% Existing nonzero fe_ppm values are preserved.
% Only rows where fe_ppm is NaN or 0 are back-filled.

if ~ismember('fe_ppm', PetroANT_DCO.Properties.VariableNames)
    PetroANT_DCO.fe_ppm = nan(height(PetroANT_DCO),1);
end

fe_ppm = PetroANT_DCO.fe_ppm;
fe_from_sum = nan(height(PetroANT_DCO),1);

has_fe2o3 = ismember('fe2o3', PetroANT_DCO.Properties.VariableNames);
has_feo   = ismember('feo',   PetroANT_DCO.Properties.VariableNames);

if has_fe2o3
    fe2o3_part = PetroANT_DCO.fe2o3 .* 1e4 .* (2 * MolMass('Fe') / MolMass('Fe2O3'));
else
    fe2o3_part = zeros(height(PetroANT_DCO),1);
end

if has_feo
    feo_part = PetroANT_DCO.feo .* 1e4 .* (1 * MolMass('Fe') / MolMass('FeO'));
else
    feo_part = zeros(height(PetroANT_DCO),1);
end

% Treat missing oxide values as zero contribution to the sum
fe2o3_part(isnan(fe2o3_part)) = 0;
feo_part(isnan(feo_part))     = 0;

fe_from_sum = fe2o3_part + feo_part;

% Only fill where there is at least some Fe oxide information
hasFeSource = false(height(PetroANT_DCO),1);
if has_fe2o3
    hasFeSource = hasFeSource | (~isnan(PetroANT_DCO.fe2o3) & PetroANT_DCO.fe2o3 ~= 0);
end
if has_feo
    hasFeSource = hasFeSource | (~isnan(PetroANT_DCO.feo) & PetroANT_DCO.feo ~= 0);
end

fillIdx = (isnan(fe_ppm) | fe_ppm == 0) & hasFeSource;
fe_ppm(fillIdx) = fe_from_sum(fillIdx);

PetroANT_DCO.fe_ppm = fe_ppm;

%% fetot_ppm = fe from fe2o3_tot
% Existing nonzero fetot_ppm values are preserved.
% feo_tot is ignored.

if ~ismember('fetot_ppm', PetroANT_DCO.Properties.VariableNames)
    PetroANT_DCO.fetot_ppm = nan(height(PetroANT_DCO),1);
end

if ismember('fe2o3_tot', PetroANT_DCO.Properties.VariableNames)
    fetot_ppm = PetroANT_DCO.fetot_ppm;
    fe_from_fe2o3tot = PetroANT_DCO.fe2o3_tot .* 1e4 .* (2 * MolMass('Fe') / MolMass('Fe2O3'));

    fillIdx = (isnan(fetot_ppm) | fetot_ppm == 0) & ...
              ~isnan(PetroANT_DCO.fe2o3_tot) & PetroANT_DCO.fe2o3_tot ~= 0;

    fetot_ppm(fillIdx) = fe_from_fe2o3tot(fillIdx);
    PetroANT_DCO.fetot_ppm = fetot_ppm;
end