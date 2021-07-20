import subprocess

print('Create training dataset for GraphProt2...')
# graphprot2 gt --in test/positives.fa --neg-in test/negatives.fa --out test_gt_out --report
process = subprocess.Popen('graphprot2', 'gt', '--in',
    f'{snakemake.input["positive"]}', '--neg-in', f'{snakemake.input["negative"]}',
    '--out', f'{snakemake.params[0]}gt_out', '--report', stdout=subprocess.PIPE,
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

print('Dataset created!')


print('Starts training GraphProt2...')
# graphprot2 train --in test_gt_out --out test_train_out
process = subprocess.Popen('graphprot2', 'train', '--in',
    f'{snakemake.params[0]}gt_out', '--out', f'{snakemake.params[0]}train_out', stdout=subprocess.PIPE,
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



print('Create prediction dataset for GraphProt2...')
# graphprot2 gp --in test/test_ws.fa --out test_gp_ws_out
process = subprocess.Popen('graphprot2', 'gp', '--in',
    f'{snakemake.input["positive"]}', '--out',
    f'{snakemake.params[0]}gp_out',
    , stdout=subprocess.PIPE,
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

print('Dataset created!')

print('Starts predicting with GraphProt2...')
# graphprot2 predict --in test_gp_ws_out --model-in test_train_out --out test_predict_ws_out --mode 1
process = subprocess.Popen('graphprot2', 'predict', '--in',
    f'{snakemake.params[0]}gp_out', '--model-in',
    f'{snakemake.params[0]}train_out',
    '--out', f'{snakemake.params[0]}predict_out', '--mode', '1', 
    stdout=subprocess.PIPE,
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
