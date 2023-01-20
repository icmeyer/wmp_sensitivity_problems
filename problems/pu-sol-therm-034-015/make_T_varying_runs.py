import os

import sens_helpers
import mpullse

def make_T_run_with_materials_script(template_dir, T):
    # Copy the template directory
    os.system('rm -r {:d}'.format(T))
    os.system('cp -r T_template {:d}'.format(T))

    cwd = os.getcwd()
    clean_dir = '{:d}/clean/'.format(T)
    materials_script = clean_dir + 'generate_materials.py'
    sens_helpers.replace_line(materials_script, 'TEMPERATURE', str(T))

    # Run materials script
    os.chdir(clean_dir)
    os.system('python generate_materials.py')
    os.chdir(cwd)

def submit_job(T):
    cwd = os.getcwd()
    os.chdir('{:d}'.format(T))
    os.system('sbatch cluster_run.sh')
    os.chdir(cwd)

if __name__=='__main__':
    Ts = mpullse.data.PU_TEMPS
    for T in Ts:
        make_T_run_with_materials_script('T_template', T)
        submit_job(T)

