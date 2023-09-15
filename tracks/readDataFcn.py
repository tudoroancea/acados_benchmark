import numpy as np
import os


def load_track(filename):
    track_file = os.path.join(os.path.dirname(__file__), filename)
    array = np.loadtxt(track_file)
    sref = array[:, 0]
    xref = array[:, 1]
    yref = array[:, 2]
    psiref = array[:, 3]
    kapparef = array[:, 4]
    return sref, xref, yref, psiref, kapparef
