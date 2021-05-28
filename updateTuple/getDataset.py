import subprocess, os, sys
from optparse import OptionParser, OptionValueError

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


cdir = os.getcwd()
print('current dir=', cdir)

import datetime

d_today = datetime.datetime.now()

d_today = str(d_today).split('.')[0]
d_today = d_today.replace(' ', '-').replace(':','')

print(d_today)



usage = "usage: python getDataset.py"
parser = OptionParser(usage)

parser.add_option("-p", "--path", default="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/Data_BcJpsiTauNu_data_2020-04-20-102531/", type="string", help="path", dest="path")

parser.add_option("-w", "--wfile", default="None", type="string", help="weight file", dest="wfile")

parser.add_option("-o", "--odir", default="None", type="string", help="odir", dest="odir")

parser.add_option("-c", "--chunk", default=3, type=int, help="chunk", dest="chunk")


(options, args) = parser.parse_args()

print('Path = ', options.path)
print('wfile = ', options.wfile)
print('output directory = ', options.odir)
print('chunk = ', options.chunk)

jobdir = cdir + '/job/' + options.odir

ensureDir(jobdir)

jlist=list(chunks(range(0, 1000), options.chunk))

#import pdb; pdb.set_trace()
#sys.exit(1)

ijob=0
for _jlist in jlist:

#    if ijob==1: continue
    
    print(ijob, _jlist)

#    import pdb; pdb.set_trace()
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
    infile = options.path

    os.system("cp job_template.sh " + jobscript)

    
    with open(jobscript) as f:
        data_lines = f.read()
        
    data_lines = data_lines.replace('INFILE', infile).replace('IDJ', '${jid}').replace('WFILE', options.wfile).replace('LOOP', ' '.join(str(x) for x in _jlist))
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)



    command = 'sbatch -p quick --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript

    os.system(command)

    ijob += 1
