function [mes, sms, oms, som] = apply_mental2(DAT)
% function pexp = apply_mental2(DAT)
% by Leonie Koban, 2024
% requires canlab core tools and Matlab
% computes pattern expression for mentalizing signature (mes),
% self-referential signature (sms), other-referential signature (oms), and
% self-vs-other signature (som) for any fmri_data object DAT (after rescaling using l2norm);
% Please cite the paper (Acil, et al., 2026, Nature Communications) if you
% use this code or the mentalizing signatures in your work.


DAT_rs = rescale(DAT, 'l2norm_images');

mes = apply_mask(DAT_rs, mes, 'pattern_expression');
sms = apply_mask(DAT_rs, sms, 'pattern_expression');
oms = apply_mask(DAT_rs, oms, 'pattern_expression');
som = apply_mask(DAT_rs, som, 'pattern_expression');

