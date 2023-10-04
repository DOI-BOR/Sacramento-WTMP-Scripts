
#from com.rma.io import DssFileManagerImpl
import os,time,sys#,shutil
from distutils.dir_util import copy_tree
from com.rma.model import Project
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))
import DMS_preprocess
reload(DMS_preprocess)

W2_models_for_input_copy = ['W2 Trinity Prescribed','W2 Lewiston Prescribed','W2 Whiskeytown Prescribed',
                            'W2 Shasta Prescribed','W2 Keswick Prescribed']

def backdate_W2_files_to_skip_compute(run_dir):
    study_dir = study_dir_from_run_dir(run_dir)
    
    current_time = time.time()
    modification_time = current_time - 3600*24*7  # Subtract 7 days (in seconds)
    
    for root, dirs, files in os.walk(os.path.join(study_dir,'cequal-w2')):
        for file in files:
            if not file.startswith('.'):
                file_path = os.path.join(root,file)
                print("Changing modified time: ",file_path)
                creation_time = os.path.getctime(file_path)
                os.utime(file_path,(creation_time,modification_time))

def study_dir_from_run_dir(run_dir):
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def annual_config_dirs_from_run_dir(run_dir,model_name,startyear_str):
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_name,model_name)  # don't know why this is two model_names deep!!
    annual_config_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,startyear_str)
    base_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,'base')
    return model_dir,annual_config_dir,base_dir

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )

    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    startyear_str = starttime_str[5:9]
    currentAlternative.addComputeMessage('Found start year for W2 simulations:'+startyear_str)
    endtime_str = rtw.getEndTimeString()
    endyear_str = endtime_str[5:9]
    if startyear_str != endyear_str:
        currentAlternative.addComputeMessage('WARNING: Start year ({0}) is different from end year ({1}); W2 simulations will likely fail.'.format(starttime_str, endtime_str))    

    run_dir = computeOptions.getRunDirectory()
    for W2_model in W2_models_for_input_copy:
        model_dir,annual_config_dir,base_dir = annual_config_dirs_from_run_dir(run_dir,W2_model,startyear_str)
        currentAlternative.addComputeMessage('model_dir: '+model_dir)
        currentAlternative.addComputeMessage('annual_config_dir: '+annual_config_dir)
        currentAlternative.addComputeMessage('base_dir: '+base_dir)
        if not os.path.exists(annual_config_dir):
            currentAlternative.addComputeMessage(W2_model+'- annual config not found; W2 may be configured incorrectly for this time window.')
        else:        
            # copy original W2 model alternative files to 'base' directory for safekeeping/later returning

            #if not os.path.exists(base_dir):
            #    os.makedir(base_dir)
            #base_files = os.listdir(base_dir)
            #if len(base_files) == 0:
            #    currentAlternative.addComputeMessage(W2_model+'- base files not found; copying from model folder')
            #    copy_tree(model_dir,base_dir)

            # remove all W2 model input files EXCEPT the .w2alt file
            for mfile in os.listdir(model_dir):
                if not mfile.endswith('.w2Alt') and not mfile.endswith('.w2Alt.bak'):
                    os.remove(os.path.join(model_dir,mfile))

            # copy over annual config input files 
            copy_tree(annual_config_dir,model_dir)
            currentAlternative.addComputeMessage('Copied W2 inputs file for '+startyear_str+' to '+W2_model+' model alternative folder')
            
            # now, W2 model alternative directory is confiured for startyear, ready for the W2 plugin to
            # work it's magic on the W2_con file and execture simulation    

	DMS_preprocess.preprocess_W2_5Res(currentAlternative, computeOptions)

	# remove me most of the time
	backdate_W2_files_to_skip_compute(run_dir)
    
    return True

