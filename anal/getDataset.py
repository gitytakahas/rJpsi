import subprocess, os, sys


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


cdir = os.getcwd()

import datetime

d_today = datetime.datetime.now()

d_today = str(d_today).split('.')[0]
d_today = d_today.replace(' ', '-').replace(':','')

print d_today



jobdir = cdir + '/job_' + d_today
outdir = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results/'
ensureDir(outdir)

if os.path.isdir(jobdir):
    ans = raw_input("Directory " + jobdir + " already exists. Delete?\n")
    
    if ans in ['Y', 'yes', 'Yes', 'y']:
        print 'deleted'
        os.system('rm -rf ' + jobdir)
    else:
        print 'quit...'
        sys.exit(1)

ensureDir(jobdir)


########################

syss = ['None']
#syss = []

for ii in range(10):
    for ud in ['up', 'down']:
        syss.append('hammer_ebe_e' + str(ii) + '_' + ud)
        
others=['puweight', 'muSFID', 'muSFReco', 'weight_ctau', 'br_BcJpsiDst', 'tauBr', 'tauReco', 'xgbsEff', 'BcPt']

for other in others:
    for ud in ['up', 'down']:
        syss.append(other + '_' + ud)

#print others

########################

#for year in ['2018']:
for year in ['2016', '2017', '2018']:
    for ijob, sys in enumerate(syss):

        jobscript = jobdir + '/job_' + year + '_' + sys + '.sh'
        
        os.system("cp job_template.sh " + jobscript)
            
        with open(jobscript) as f:
            data_lines = f.read()
        
        data_lines = data_lines.replace('SYSTEMATIC',sys).replace('OUTDIRECTORY', outdir).replace('YEARTOBEFILLED', year)
        
        with open(jobscript, mode="w") as f:
            f.write(data_lines)

        command = 'sbatch -p short --account=t3 --error=' + jobdir + '/err_' + year + '_' + sys + '.' + str(ijob) + ' --output=' + jobdir + '/out_' + year + '_' + sys + '.' + str(ijob) + ' ' + jobscript

        print year, sys
#        print(command)
        os.system(command)
