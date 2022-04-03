import numpy as np
from scipy.spatial.transform import Rotation

def computeCenterOfMass(indices, top):
    tM = 0.0
    com = np.array([0.0, 0.0, 0.0])

    masses = {}

    for tinfo in top.propertiesLoaded['TYPES']:
        t = tinfo[0].split()[0]
        m = float(tinfo[0].split()[1])
        masses[t] = m
        
    for i in indices:
        pos = np.asarray([float(p) for p in top.coordLoaded[i][1].split()])
        tp = top.propertiesLoaded['STRUCTURE'][i][1]
        tM = tM + masses[tp]
        com = com + pos * masses[tp]
    
    return com / tM


def rotateStructure(axis, angle, top):
    rot = Rotation.from_rotvec(angle * axis)
    for i, c in enumerate(top.coordLoaded):
        pos = np.asarray([float(p) for p in c[1].split()])
        pos = rot.apply(pos)
        pNew = [c[0], '{:.4f} {:.4f} {:.4f}'.format(pos[0], pos[1], pos[2])]
        top.coordLoaded[i] = pNew.copy()

