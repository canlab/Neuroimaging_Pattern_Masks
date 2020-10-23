% Given an fmri_data object, will return pattern expressions for the 5
% patterns in MPA2 multiaversive

function pexp_tabl = apply_multiaversive_mpa2_patterns(dat)

    pats = load_image_set('mpa2');
    
    pexps = apply_mask(dat, pats, 'pattern_expression');
    
    pexp_tabl = array2table(pexps, 'VariableNames', {'general' 'mechanical' 'thermal' 'sound' 'visual'});
    
end
    
    

