% Convert elemental *_ppm columns to *_mol columns in table T
% Assumes ppm = ug/g, so:
%   mol/g = ppm * 1e-6 / MolMass(element)
%
% Resulting columns are named like:
%   si_ppm  -> si_mol
%   fetot_ppm -> fetot_mol
%
% Existing *_mol columns are overwritten.

ppmMap = {
    'f_ppm',      'F'
    'cl_ppm',     'Cl'
    'br_ppm',     'Br'
    'i_ppm',      'I'
    'h_ppm',      'H'
    'c_ppm',      'C'
    'n_ppm',      'N'
    'p_ppm',      'P'
    's_ppm',      'S'
    'al_ppm',     'Al'
    'as_ppm',     'As'
    'ag_ppm',     'Ag'
    'au_ppm',     'Au'
    'b_ppm',      'B'
    'ba_ppm',     'Ba'
    'be_ppm',     'Be'
    'bi_ppm',     'Bi'
    'ca_ppm',     'Ca'
    'cd_ppm',     'Cd'
    'ce_ppm',     'Ce'
    'co_ppm',     'Co'
    'cr_ppm',     'Cr'
    'cs_ppm',     'Cs'
    'cu_ppm',     'Cu'
    'dy_ppm',     'Dy'
    'er_ppm',     'Er'
    'eu_ppm',     'Eu'
    'fe_ppm',     'Fe'
    'ga_ppm',     'Ga'
    'gd_ppm',     'Gd'
    'ge_ppm',     'Ge'
    'hf_ppm',     'Hf'
    'hg_ppm',     'Hg'
    'ho_ppm',     'Ho'
    'in_ppm',     'In'
    'ir_ppm',     'Ir'
    'k_ppm',      'K'
    'la_ppm',     'La'
    'li_ppm',     'Li'
    'lu_ppm',     'Lu'
    'mg_ppm',     'Mg'
    'mn_ppm',     'Mn'
    'mo_ppm',     'Mo'
    'na_ppm',     'Na'
    'nd_ppm',     'Nd'
    'ni_ppm',     'Ni'
    'nb_ppm',     'Nb'
    'os_ppm',     'Os'
    'pa_ppm',     'Pa'
    'pb_ppm',     'Pb'
    'pd_ppm',     'Pd'
    'pm_ppm',     'Pm'
    'pr_ppm',     'Pr'
    'pt_ppm',     'Pt'
    'rb_ppm',     'Rb'
    're_ppm',     'Re'
    'rh_ppm',     'Rh'
    'ru_ppm',     'Ru'
    'sb_ppm',     'Sb'
    'sc_ppm',     'Sc'
    'se_ppm',     'Se'
    'si_ppm',     'Si'
    'sm_ppm',     'Sm'
    'sn_ppm',     'Sn'
    'sr_ppm',     'Sr'
    'ta_ppm',     'Ta'
    'tb_ppm',     'Tb'
    'te_ppm',     'Te'
    'th_ppm',     'Th'
    'ti_ppm',     'Ti'
    'tl_ppm',     'Tl'
    'tm_ppm',     'Tm'
    'w_ppm',      'W'
    'v_ppm',      'V'
    'u_ppm',      'U'
    'y_ppm',      'Y'
    'yb_ppm',     'Yb'
    'zn_ppm',     'Zn'
    'zr_ppm',     'Zr'
    'fetot_ppm',  'Fe'
};

for i = 1:size(ppmMap,1)
    ppmVar = ppmMap{i,1};
    elem   = ppmMap{i,2};

    if ~ismember(ppmVar, PetroANT_DCO.Properties.VariableNames)
        warning('Missing column %s',ppmVar);
        continue
    end

    molVar = regexprep(ppmVar, '_ppm$', '_mol');
    PetroANT_DCO.(molVar) = PetroANT_DCO.(ppmVar) .* 1e-6 ./ MolMass(elem);
end