# Regenerate the circular brain thumbnails embedded in the two overview figures.
# Sources are the per-folder png_images/*_surface.png (and *_isosurface.png)
# renders already checked into the repo. Output -> docs/figure_source/_thumbs/.
#
#   python3 build_thumbnails.py
#
# Requires Pillow (pip install pillow). The _thumbs/ dir is regenerable and is
# git-ignored; the finished SVGs embed these images as base64.
import os
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # repo root
OUT_OV = os.path.join(HERE, "_thumbs", "overview")
OUT_SIG = os.path.join(HERE, "_thumbs", "signatures")
os.makedirs(OUT_OV, exist_ok=True)
os.makedirs(OUT_SIG, exist_ok=True)

# surface renders are a 2x2 montage of brain views; this crops the top-left
# lateral view. isosurface renders are a single brain -> no crop (quad=None).
QUAD = (40, 25, 378, 300)


def whiten(im, thr=238):
    im = im.convert("RGBA")
    px = im.load()
    w, h = im.size
    for y in range(h):
        for x in range(w):
            r, g, b, a = px[x, y]
            if r >= thr and g >= thr and b >= thr:
                px[x, y] = (r, g, b, 0)      # white background -> transparent
    return im


def process(src, out, quad=QUAD, size=200, pad=0.10):
    im = Image.open(os.path.join(ROOT, src)).convert("RGBA")
    if quad:
        im = im.crop(quad)
    im = whiten(im)
    bb = im.getbbox()
    if bb:
        im = im.crop(bb)
    w, h = im.size
    side = int(max(w, h) * (1 + pad))
    canvas = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    canvas.paste(im, ((side - w) // 2, (side - h) // 2), im)
    canvas.resize((size, size), Image.LANCZOS).save(out)


# ---- five thumbnails for the repository-overview figure ----
OVERVIEW = {
    "sig":           ("Multivariate_signature_patterns/2017_Woo_SIIPS1/png_images/SIIPS1_weighted_mean_surface.png", QUAD),
    "atlas":         ("Atlases_and_parcellations/2024_CANLab_atlas/png_images/canlab2024_fine_fmriprep20_2mm_isosurface.png", None),
    "meta_kober":    ("CANlab_Meta_analysis_maps/2008_Kober_Emotion_163_studies/png_images/Kober2008_EmoActivation_FWE_all_surface.png", QUAD),
    "neurosynth_iso":("Neurosynth_maps/mkda/png_images/Neurosynth_MKDA_Activation_proportion_isosurface.png", None),
    "study":         ("Individual_study_maps/2024_Bo_EmotionRegulation_BayesFactor/png_images/Bo2024_CommonAppraisal_surface.png", QUAD),
}

# ---- one thumbnail per signature chip in the taxonomy figure ----
# key -> repo-relative *_surface.png. Signatures with no surface render in the
# repo (somatovisceral, Matthewson SCR, saCPM) are drawn with a vector glyph in
# build_taxonomy.py and are intentionally absent here.
SIG = "Multivariate_signature_patterns"
SIGNATURES = {
    "nps":        f"{SIG}/2017_Lopez_Sola_Fibromyalgia/png_images/LopezSola2017_NPSp_surface.png",
    "siips":      f"{SIG}/2017_Woo_SIIPS1/png_images/SIIPS1_weighted_mean_surface.png",
    "feps":       f"{SIG}/2024_FEPS_Facial_Expressions_of_Pain_Signature/png_images/FEPS_mean_xval_weights_surface.png",
    "cpdm":       f"{SIG}/2020_Geuter_pain_multivariate_mediation_PDM/png_images/Geuter2020_cPDM_combined_surface.png",
    "vps":        f"{SIG}/2016_Krishnan_eLife_VPS/png_images/VPS_unthresholded_surface.png",
    "genvic":     f"{SIG}/2020_Zhou_general_vicarious_pain/png_images/Zhou2020_General_unthresh_surface.png",
    "placebo":    f"{SIG}/2011_Wager_JNeuro_placebo_prediction/png_images/PlaceboPredict_Anticipation_surface.png",
    "paincog":    f"{SIG}/2020_Silvestrini_Rainville_Pain_CogControl_interaction_aMCC/png_images/Silvestrini2020_dACC_pain_surface.png",
    "painval":    f"{SIG}/2022_coll_pain_monetary_reward_decision_value/png_images/Coll2022_painvalue_unthresh_surface.png",
    "gsr":        f"{SIG}/2016_Eisenbarth_JNeuro_autonomic_patterns/png_images/Eisenbarth2016_GSR_unthresh_surface.png",
    "hr":         f"{SIG}/2016_Eisenbarth_JNeuro_autonomic_patterns/png_images/Eisenbarth2016_HR_unthresh_surface.png",
    "threat":     f"{SIG}/2018_Reddan_Threat_Conditioning_ImEx/png_images/Reddan2018_CSplus_unthresh_surface.png",
    "vifs":       f"{SIG}/2021_Zhou_Subjective_Fear/png_images/Zhou2021_VIFS_surface.png",
    "mpa2":       f"{SIG}/2021_Ceko_MPA2_multiaversive/png_images/MPA2_General_surface.png",
    "pines":      f"{SIG}/2015_Chang_PLoSBiology_PINES/png_images/PINES_Rating_Weights_LOSO_surface.png",
    "bpls_neg":   f"{SIG}/2015_Kragel_emotionClassificationBPLS/png_images/Kragel2015_BPLS_fearful_surface.png",
    "schema_neg": f"{SIG}/2019_Kragel_Emotion_Schemas/png_images/Kragel2019_Fear_surface.png",
    "distress":   f"{SIG}/2017_Ashar_care_distress/png_images/Ashar2017_distress_unthresh_surface.png",
    "guilt":      f"{SIG}/2019_Yu_Koban_Guilt/png_images/Yu2020_Guilt_SVM_surface.png",
    "rejection":  f"{SIG}/2015_Woo_NatureComms_Rejection/png_images/Rejection_dpsp_weights_surface.png",
    "brs":        f"{SIG}/2023_Speer_Brain_Reward_Signature_BRS/png_images/Speer2023_BRS_surface.png",
    "ncs":        f"{SIG}/2022_Koban_NCS_Craving/png_images/Koban2023_NCS_general_surface.png",
    "basic":      f"{SIG}/2021_vantHoff_BASIC_sexual_image_classifier/png_images/vantHoff2021_BASIC_surface.png",
    "rewval":     f"{SIG}/2022_coll_pain_monetary_reward_decision_value/png_images/Coll2022_moneyvalue_unthresh_surface.png",
    "care":       f"{SIG}/2017_Ashar_care_distress/png_images/Ashar2017_care_unthresh_surface.png",
    "bpls_pos":   f"{SIG}/2015_Kragel_emotionClassificationBPLS/png_images/Kragel2015_BPLS_amused_surface.png",
    "schema_pos": f"{SIG}/2019_Kragel_Emotion_Schemas/png_images/Kragel2019_Joy_surface.png",
    "stroop":     f"{SIG}/2020_Silvestrini_Rainville_Pain_CogControl_interaction_aMCC/png_images/Silvestrini2020_Stroop_WB_surface.png",
    "mfc":        f"{SIG}/2018_Kragel_MFC_Generalizability/png_images/Kragel2018_bPLS_Wholebrain_Cognitive_Control_surface.png",
    "mentalizing":f"{SIG}/2026_Acil_Mentalizing_Self_Other/png_images/Acil2026_MS_unthresh_surface.png",
    "fibro":      f"{SIG}/2017_Lopez_Sola_Fibromyalgia/png_images/LopezSola2017_FM_pain_surface.png",
    "backpain":   f"{SIG}/2019_Lee_JPain_backpain/png_images/Lee2019_S1_pairedSVM_surface.png",
    "pifonem":    f"{SIG}/2026_Murillo_PiFoneM/png_images/Murillo2026_PiFoneM_unthresh_surface.png",
}


def main():
    for key, (src, quad) in OVERVIEW.items():
        process(src, os.path.join(OUT_OV, key + ".png"), quad=quad, size=240, pad=0.12)
    print(f"overview thumbnails: {len(OVERVIEW)} -> {OUT_OV}")
    missing = []
    for key, src in SIGNATURES.items():
        if not os.path.exists(os.path.join(ROOT, src)):
            missing.append(key)
            continue
        # signature chips are small; a surface render is always a 2x2 montage
        quad = None if "isosurface" in src else QUAD
        process(src, os.path.join(OUT_SIG, key + ".png"), quad=quad, size=168, pad=0.12)
    print(f"signature thumbnails: {len(SIGNATURES) - len(missing)} -> {OUT_SIG}")
    if missing:
        print("MISSING sources for:", missing)


if __name__ == "__main__":
    main()
