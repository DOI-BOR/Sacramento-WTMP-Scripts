#version 2.0
#modified 03-28-2023 by Scott Burdick-Yahya

from com.rma.model import Project
import os

def fixInputLocationFpart(currentAlternative, tspath):
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    tspath = tspath.split('/')
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

def appendAPart(current_path, ApartAppend):
    tspath = tspath.split('/')
    Apart = tspath[1]
    if len(Apart) == 0:
        new_Apart = ApartAppend
    else:
        new_Apart = Apart + '_' + ApartAppend
    tspath[1] = new_Apart
    tspath = '/'.join(tspath)
    return tspath

def getDataLocationDSSInfo(location, currentAlternative, computeOptions):
    if location.isLinkedToPreviousModel():
        tspath = str(currentAlternative.loadTimeSeries(location))
        tspath = fixInputLocationFpart(currentAlternative, tspath)
        dsspath = computeOptions.getDssFilename()
    else:
        tspath = location.getLinkedToLocation().getDssPath()
        rundir = Project.getCurrentProject().getProjectDirectory()
        dsspath = location.getLinkedToLocation().get_dssFile()
        dsspath = os.path.join(rundir, dsspath)
    return tspath, dsspath
