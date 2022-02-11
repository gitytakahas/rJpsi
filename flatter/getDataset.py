import subprocess, os, sys
from optparse import OptionParser, OptionValueError

prefix = 'gsiftp://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/'


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


cdir = os.getcwd()


import datetime

d_today = datetime.datetime.now()

d_today = str(d_today).split('.')[0]
d_today = d_today.replace(' ', '-').replace(':','')

print d_today



usage = "usage: python getDataset.py"
parser = OptionParser(usage)

parser.add_option("-f", "--file", default="RJpsi_ParkingBPH1_2020-01-08-180946_20200108180945_dz0_fsig3_vprob0p1/ParkingBPH1/ParkingBPH1_20200108180945_dz0_fsig3_vprob0p1/200108_171159/", type="string", help="file", dest="file")
parser.add_option("-c", "--chunk", default=15, type="int", help="chunk", dest="chunk")
parser.add_option("-a", "--analysis", default="BsTauTau", type="string", help="analysis channel", dest="analysis")
parser.add_option("-t", "--type", default="data", type="string", help="type", dest="type")
parser.add_option("-n", "--name", default="None", type="string", help="name", dest="name")
parser.add_option("-s", "--select", default="None", type="string", help="select", dest="select")
parser.add_option("-e", "--exclude", default="None", type="string", help="exclude", dest="exclude")
parser.add_option("-y", "--year", default="UL2017", type="string", help="year", dest="year")
parser.add_option("-o", "--odir", default="job", type="string", help="output dir", dest="odir")
parser.add_option("-j", "--jdir", default="job", type="string", help="job output dir", dest="jdir")
parser.add_option("-p", "--priority", default="pt", type="string", help="priority", dest="priority")

(options, args) = parser.parse_args()

print 'File = ', options.file
print 'Chunks = ', options.chunk




#RJpsi_ParkingBPH1_2020-01-19-130139_20200119130138_BsPhiTauTau_mutau/ParkingBPH1/ParkingBPH1_Run2018A-05May2019-v1_20200119130138_BsPhiTauTau_mutau/200119_120740/0000

out = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.file], 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT)


listoffiles = []



for line in iter(out.stdout.readline,''):

   line = line.rstrip()
   

   if len( line.split())!=9: continue

   ndir = line.rstrip().split()[8]

   print 'Adding ...', prefix + '/' + options.file + '/' + ndir

   out2 = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.file + '/' + ndir],
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.STDOUT)

   for line2 in iter(out2.stdout.readline,''):

      line2 = line2.rstrip()

#      print line2

      if len( line2.split())!=9: continue


      nndir = line2.split()[8]

      
      out3 = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.file + '/' + ndir + '/' + nndir],
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.STDOUT)
      
      

      for line3 in iter(out3.stdout.readline,''):
          
          line3 = line3.rstrip()
          
          if len( line3.split())!=9: continue

          nnndir = line3.split()[8]

          out4 = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.file + '/' + ndir + '/' + nndir + '/' + nnndir],
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT)

          for line4 in iter(out4.stdout.readline,''):

              line4 = line4.rstrip()

              if len( line4.split())!=9: continue


              nnnndir = line4.split()[8]
              
      
              out5 = subprocess.Popen(['uberftp', '-ls', prefix + '/' + options.file + '/' + ndir + '/' + nndir + '/' + nnndir + '/' + nnnndir],
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.STDOUT)
              
              for line5 in iter(out5.stdout.readline,''):

                  line5 = line5.rstrip()
                  
                  if len( line5.split())!=9: continue

                  if line5.find('.root')==-1: continue

                  file = line5.split()[8]

                  filename = prefix.replace('gsiftp', 'root') + '/' + options.file + '/' + ndir + '/' + nndir + '/' + nnndir + '/' + nnnndir + '/' + file
                  print '\t', filename 
                  listoffiles.append(filename)



#sys.exit(1)

print len(listoffiles), 'files detected'

if options.select!="None":
    listoffiles = [f for f in listoffiles if f.find(options.select)!=-1]
    print len(listoffiles), 'selected with string, ', options.select

if options.exclude!="None":
    listoffiles = [f for f in listoffiles if f.find(options.exclude)==-1]
    print len(listoffiles), 'selected with excluded string, ', options.exclude


print '-'*80
for file_ in listoffiles:
    print '->', file_
print '-'*80

#print len(listoffiles), 'files are detected'

listoffiles = list(chunks(listoffiles,options.chunk))

print len(listoffiles), 'chunks are created'

#ans = raw_input("Is this OK?\n")
    
#if ans not in ['Y', 'yes', 'Yes']:
#    sys.exit(1)


jobdir = cdir + '/' + options.jdir + '/' + options.name #+ '_' + options.analysis + '_' + options.type
outdir = options.odir + '/' + options.jdir + '/' + options.name
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


for ijob, files in enumerate(listoffiles):
#    print ijob

#    if ijob==2: break
    input = ','.join(files)

    # Now, write the shell script into job directory
    
    jobscript = jobdir + '/job_' + str(ijob) + '.sh'
    outfile = outdir + '/Myroot_' + str(ijob) + '.root'

    os.system("cp job_template.sh " + jobscript)

    
    with open(jobscript) as f:
        data_lines = f.read()
        
    data_lines = data_lines.replace('PROCESS',options.analysis).replace('INFILE', input).replace('OUTFILE', outfile).replace('TYPE', options.type).replace('YEAR', options.year).replace('IDJ', str(ijob)).replace('PRIORITY', options.priority).replace('ONAME', options.name)
        
    with open(jobscript, mode="w") as f:
        f.write(data_lines)

    command = 'sbatch -p standard --account=t3 --error=' + jobdir + '/err.' + str(ijob) + ' --output=' + jobdir + '/out.' + str(ijob) + ' ' + jobscript
    print(command)
    os.system(command)
