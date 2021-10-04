import random
import os
import argparse


def main(input, output):
    random.seed(42)

    tmp_dir = 'shuffle/'
    os.mkdir(tmp_dir)

    n = 10 # broj tmp fajlova koji se prave, menjati po potrebi

    paths = [f'{tmp_dir}{i}' for i in range(n)]
    files = [open(path, 'w') for path in paths]

    with open(input, 'r') as input_file:
        for row in input_file:
            current = random.randint(0, n - 1)
            files[current].write(row)
            if row[-1] != '\n':
                files[current].write('\n')

    for file in files:
        file.close()

    with open(output, 'w') as output_file:
        for path in paths:
            with open(path, 'r') as file:
                lines = file.readlines()
                random.shuffle(lines)
                output_file.writelines(lines)

    for path in paths:
        os.remove(path)
    os.rmdir(tmp_dir)

main('data/interim/novi.csv', 'data/interim/testing.csv')

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(
#         description='Shuffle file lines.'
#     )
#
#     parser.add_argument(
#         '--input', '-i',
#         dest='input',
#         action='store',
#         type=str,
#         required=True,
#         help='Path to file to shuffle.'
#     )
#     parser.add_argument(
#         '--output', '-o',
#         dest='output',
#         action='store',
#         type=str,
#         required=True,
#         help='Path to file to save the shuffled lines.'
#     )
#
#     args = parser.parse_args()
#     main(args)
