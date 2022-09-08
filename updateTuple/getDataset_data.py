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

parser.add_option("-p", "--path", default="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_multiple_v2/Data_2018/", type="string", help="path", dest="path")

parser.add_option("-o", "--odir", default="None", type="string", help="odir", dest="odir")

parser.add_option("-c", "--chunk", default=10, type=int, help="chunk", dest="chunk")

parser.add_option("-n", "--name", default="None", type="string", help="name", dest="name")

parser.add_option("-m", "--model", default="None", type="string", help="model", dest="model")

parser.add_option("-j", "--jdir", default="None", type="string", help="jdir", dest="jdir")

parser.add_option('-s', '--signal', action="store_true", default=False, dest='signal')

(options, args) = parser.parse_args()

print('Path = ', options.path)
print('output directory = ', options.odir)
print('job output directory = ', options.jdir)
print('chunk = ', options.chunk)
print('name = ', options.name)

jobdir = cdir + '/job/' + options.jdir

ensureDir(jobdir)

#jlist=list(chunks(range(-1, 1000), options.chunk))

#import pdb; pdb.set_trace()
#sys.exit(1)


from os import listdir
from os.path import isfile, join
jlist = [f for f in listdir(options.path) if isfile(join(options.path, f)) and f.find('.root')!=-1 and f.find('Myroot')!=-1 and f!="Myroot.root" and f!="Myroot_data_2018.root" and f.find('training')==-1 and f.find('analysis')==-1 and f.find('data')==-1]

if options.signal:
    jlist = [f for f in listdir(options.path) if isfile(join(options.path, f)) and f.find('hammer')!=-1]


#import pdb; pdb.set_trace()

print(jlist)
#print len(onlyfiles), 'files detected'

#jlist=list(chunks(onlyfiles, options.chunk))

print(len(jlist), 'jobs created')

#sys.exit(1)

ijob=0
for _jlist in jlist:
    
    print(ijob, _jlist)

#    if ijob == 1: break

#    import pdb; pdb.set_trace()
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
#    infile = options.path

    os.system("cp job_template_data.sh " + jobscript)

    
    with open(jobscript) as f:
        data_lines = f.read()
        
#    data_lines = data_lines.replace('LOOP', ' '.join(options.path + '/' + str(x) for x in _jlist)).replace('ODIR', options.odir).replace('PNAME', options.name).replace('OMODEL', options.model).replace('IDJ', str(ijob))
    data_lines = data_lines.replace('INFILE', options.path + '/' + _jlist).replace('ODIR', options.odir).replace('PNAME', options.name + '_' + str(ijob)).replace('OMODEL', options.model).replace('OUTFILE','tmp_' + str(ijob) + '.root')
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)



    command = 'sbatch -p short --mem=4G --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript

    os.system(command)

    ijob += 1


print(ijob, 'created for', len(jlist), 'files')
