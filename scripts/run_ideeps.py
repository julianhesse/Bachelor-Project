import os
import subprocess

os.chdir('./methods/iDeepS')
if not os.path.isdir('./EDeN'):
    subprocess.run(['unzip', 'EDeN.zip'])
os.chdir('./EDeN')

# python setup.py build_ext -I$CONDA_PREFIX/include/openbabel3 -L$CONDA_PREFIX/lib
subprocess.run(['python', 'setup.py', 'build_ext', '-I$CONDA_PREFIX/include/openbabel3', '-L$CONDA_PREFIX/lib'])
subprocess.run(['python', 'setup.py', 'install'])
os.chdir('../')

def clean_up():
    print('Clean up ...')
    if os.path.isfile('structure.gz'):
        os.remove('structure.gz')
        print('Cleaned up structures.gz')
    else:
        print('Nothing to clean up')

    data_file_dir = os.path.dirname('../../{snakemake.input["train"]}')
    if os.path.isfile(data_file_dir + '/structure.gz'):
        os.remove(data_file_dir + '/structure.gz')
        print('Cleaned up %s/structures.gz' % data_file_dir)
    else:
        print('Nothing to clean up elsewhere')

# clean up
clean_up()

print('Starts training iDeepS...')
# python ideeps.py --train=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/training_sample_0/sequences.fa.gz --model_dir=models
process = subprocess.Popen(['python', 'ideeps.py', '--train=True',
    f'--data_file=../../{snakemake.input["train"]}', f'--model_dir=../../{snakemake.params[0]}'], stdout=subprocess.PIPE,
    universal_newlines=True)

while True:
    output = process.stdout.readline()
    print(output.strip())
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        for output in process.stdout.readlines():
            print(output.strip())
        break

print('\nTraining finished!\n')

# clean up again
clean_up()

data_file_dir = os.path.dirname('../../{snakemake.input["train"]}')
if os.path.isfile(data_file_dir + '/structure.gz'):
    os.remove(data_file_dir + '/structure.gz')
    print('Cleaned up %s/structures.gz' % data_file_dir)
else:
    print('Nothing to clean up')

# python ideeps.py --predict=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --out_file=YOUR_OUTFILE
print('Starts predicting with iDeepS...')
process = subprocess.Popen(['python', 'ideeps.py', '--predict=True',
    f'--data_file=../../{snakemake.input["test"]}', f'--model_dir=../../{snakemake.params[0]}',
    f'--out_file=../../{snakemake.output["prediction"]}'], stdout=subprocess.PIPE,
    universal_newlines=True)

while True:
    output = process.stdout.readline()
    print(output.strip())
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        for output in process.stdout.readlines():
            print(output.strip())
        break

print('Prediction finished!')

# clean up again
clean_up()
