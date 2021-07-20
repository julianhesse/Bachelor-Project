import subprocess

print('Starts training iDeepS...')
# python ideeps.py --train=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/training_sample_0/sequences.fa.gz --model_dir=models
process = subprocess.Popen('python', 'ideeps.py', '--train=True',
    f'--data_file={snakemake.input[0]}', f'--model_dir={snakemake.params[0]}', stdout=subprocess.PIPE,
    universal_newlines=True)

while True:
    output = process.stdout.readline()
    print(output.strip())
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        for output in precess.stdout.readlines():
            print(output.strip())
        break

print('Training finished!')


# python ideeps.py --predict=True --data_file=datasets/clip/10_PARCLIP_ELAVL1A_hg19/30000/test_sample_0/sequences.fa.gz --model_dir=models --out_file=YOUR_OUTFILE
print('Starts predicting with iDeepS...')
process = subprocess.Popen('python', 'ideeps.py', '--predict=True',
    f'--data_file={snakemake.input[0]}', f'--model_dir={snakemake.params[0]}', stdout=subprocess.PIPE,
    universal_newlines=True)

while True:
    output = process.stdout.readline()
    print(output.strip())
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # read rest of output
        for output in precess.stdout.readlines():
            print(output.strip())
        break

print('Prediction finished!')
