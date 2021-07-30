import subprocess
import os


## Test if GraphProt2 image needs to be created
print('\nChecking if image for GraphProt2 exists...')
if os.path.isdir('/var/tmp/graphprot2'):
    print('Image already build!')
else:
    os.chdir('./methods/GraphProt2')
    print('Building the image!')

    # Build image with:
    # ch-build -t graphprot2 .
    process = subprocess.Popen(['ch-build', '-t', 'graphprot2', '.'],
        stdout=subprocess.PIPE, universal_newlines=True)

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

    print('\nPacking build into tar')
    # Create tar:
    # ch-builder2tar graphprot2 /var/tmp
    process = subprocess.Popen(['ch-builder2tar', 'graphprot2', '/var/tmp'],
        stdout=subprocess.PIPE, universal_newlines=True)

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

    print('\nUnpacking tar to destination')
    # Unpack tar:
    # ch-tar2dir /var/tmp/graphprot2.tar.gz /var/tmp
    process = subprocess.Popen(['ch-tar2dir', '/var/tmp/graphprot2.tar.gz', '/var/tmp'],
        stdout=subprocess.PIPE, universal_newlines=True)

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

    print('\nImage is created successfully!')
    os.chdir('../../')

print('Running GraphProt2...')
# Run graphprot2:
# ch-run -b path_to_project:/data /var/tmp/graphprot2 \
# -- ./methods/Graphprot2/wrapper.sh path_to_positive_fasta \
# path_to_negative_fasta path_to_test_fasta path_to_ouput_folder
positive_file = snakemake.input['positive']
negative_file = snakemake.input['negative']
test_file = snakemake.input['test']
out_folder = snakemake.params[0]
process = subprocess.Popen(['ch-run', '-b', f'{os.getcwd()}:/data', '/var/tmp/graphprot2',
    './methods/Graphprot2/wrapper.sh', positive_file, negative_file,
    test_file, out_folder],
    stdout=subprocess.PIPE, universal_newlines=True)

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

print('\nGraphProt2 did run successfully!')

""" Old way of execution
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
        for output in process.stdout.readlines():
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
        for output in process.stdout.readlines():
            print(output.strip())
        break

print('Training finished!')



print('Create prediction dataset for GraphProt2...')
# graphprot2 gp --in test/test_ws.fa --out test_gp_ws_out
process = subprocess.Popen('graphprot2', 'gp', '--in',
    f'{snakemake.input["positive"]}', '--out',
    f'{snakemake.params[0]}gp_out',
    stdout=subprocess.PIPE, universal_newlines=True)

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
        for output in process.stdout.readlines():
            print(output.strip())
        break

print('Prediction finished!')
"""
