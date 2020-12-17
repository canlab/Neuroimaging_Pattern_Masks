# -*-coding:Utf-8 -*

import subprocess


def ctrl(opt_name, opt_value):
    if opt_value == None:
        return ""
    else:
        return " " + opt_name + " " + str(opt_value)

    # We can use the mask of the brain resulting from BrainExtraction
    # in Atropos to define the regions we want to segment.


def segment(dim, anatomical, mask=None, nb_class=None,
            out_pref=None, priors=None):
    opt = "-d " + str(dim) \
        + " -a " + anatomical \
        + ctrl("-x", mask) \
        + ctrl("-c", nb_class) \
        + ctrl("-o", out_pref) \
        + ctrl("-p", priors)
    s = subprocess.Popen(["antsAtroposN4.sh", opt], stdout=subprocess.PIPE)
    out, err = s.communicate()
    if err != None:
        print("Error in antsAtroposN4.sh : \n" + err)
    else:
        print(out)
