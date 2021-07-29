import os
import subprocess

os.chdir('./methods/iDeepS/EDeN')

# python setup.py build_ext -I$CONDA_PREFIX/include/openbabel3 -L$CONDA_PREFIX/lib
subprocess.run(['python', 'setup.py', 'build_ext', '-I$CONDA_PREFIX/include/openbabel3', '-L$CONDA_PREFIX/lib'])
subprocess.run(['python', 'setup.py', 'install'])
os.chdir('../')

print('Starts training iDeepS...')
# python ideeps.py --train=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/training_sample_0/sequences.fa.gz --model_dir=models
process = subprocess.Popen(['python', 'ideeps.py', '--train=True',
    f'--data_file=../../{snakemake.input[0]}', f'--model_dir=../../{snakemake.params[0]}'], stdout=subprocess.PIPE,
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

print('Training finished!')


# python ideeps.py --predict=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --out_file=YOUR_OUTFILE
print('Starts predicting with iDeepS...')
process = subprocess.Popen(['python', 'ideeps.py', '--predict=True',
    f'--data_file=../../{snakemake.input[0]}', f'--model_dir=../../{snakemake.params[0]}',
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

# clean up
print('Clean up ...')
subprocess.run(['rm', 'structure.gz'])
print('Cleaned up structures.gz')
