import os
import subprocess
import sys

input_file = sys.argv[1]
test_file = sys.argv[2]
prediction_file = sys.argv[3]
out_dir = sys.argv[4]
fold = sys.argv[5]


os.chdir('./methods/ideeps_os')
if not os.path.isdir('./EDeN'):
    # subprocess.run(['unzip', 'EDeN.zip'])
    print subprocess.Popen(['unzip', 'EDeN.zip'], stdout=subprocess.PIPE).stdout.read()
os.chdir('./EDeN')

# python setup.py build_ext -I$CONDA_PREFIX/include/openbabel3 -L$CONDA_PREFIX/lib
subprocess.Popen(['python', 'setup.py', 'build_ext', '-I$CONDA_PREFIX/include/openbabel3', '-L$CONDA_PREFIX/lib'])
subprocess.Popen(['python', 'setup.py', 'install'])
os.chdir('../')

def clean_up():
    print('Clean up ...')
    if os.path.isfile('structure.gz'):
        os.remove('structure.gz')
        print('Cleaned up structures.gz')
    else:
        print('Nothing to clean up')

    model_dir = '../../' + out_dir
    if os.path.isfile(model_dir + 'structure.gz'):
        os.remove(model_dir + 'structure.gz')
        print('Cleaned up %sstructures.gz' % model_dir)
    else:
        print('Nothing to clean up at %s' % model_dir)

# clean up
clean_up()

# create theano compile folder and necessary env
theano_folder = "%s/../../.theano-%s" % (os.getcwd(), fold)
if not os.path.isdir(theano_folder):
    subprocess.Popen(['mkdir', theano_folder])
my_env = os.environ.copy()
my_env['THEANO_FLAGS'] = "base_compiledir=%s" % theano_folder

print '\nStarts training iDeepS...'
# python ideeps.py --train=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/training_sample_0/sequences.fa.gz --model_dir=models
process = subprocess.Popen(['python', 'ideeps.py', '--train=True',
    '--data_file=../../%s' % input_file, '--model_dir=../../%s' % out_dir], stdout=subprocess.PIPE,
    universal_newlines=True, env=my_env)

while True:
    output = process.stdout.readline()
    out = output.strip()
    if out != '':
        print out
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        # for output in process.stdout.readlines():
            # print(output.strip())
        break

print '\nTraining finished!\n'

# clean up again
clean_up()

# python ideeps.py --predict=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --out_file=YOUR_OUTFILE
print '\nStarts predicting with iDeepS...'
process = subprocess.Popen(['python', 'ideeps.py', '--predict=True',
    '--data_file=../../%s' % test_file, '--model_dir=../../%s' % out_dir,
    '--out_file=../../%s' % prediction_file], stdout=subprocess.PIPE,
    universal_newlines=True, env=my_env)

while True:
    output = process.stdout.readline()
    out = output.strip()
    if out != '':
        print out
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        # for output in process.stdout.readlines():
            # print(output.strip())
        break

print '\nPrediction finished!\n'

# clean up again
clean_up()
