import os
import json
import fnmatch


def search_fname(root, pattern):

    file_paths = []
    file_list = os.listdir(root)
    for filename in file_list:
        if fnmatch.fnmatch(filename, pattern):
            file_paths.append(os.path.join(root, filename))
    if not file_paths: raise FileNotFoundError("Matching files not found.")
    else: print("Found matching files:", file_paths )

    control=0
    chosen_file_path = file_paths
    opt_read = input('Do you want to read params from file? [[y]/n]')
    if opt_read=='y' or opt_read=='':
        if len(file_paths)==1:
            print(f'You have chosen scapolottina: {file_paths}')
            chosen_file_path = file_paths[0]
        else:
            opt_params_year = input('Which year do you want to use for static parameters A,B,D?')
            print('Chosen file of year:', opt_params_year)
            pattern = '*_'+opt_params_year+'_*'
            for filename in file_paths:
                if fnmatch.fnmatch(filename, pattern):
                    if control: raise FileExistsError('Multiple files found: please check for duplicates.')
                    else:
                        print(f'You have chosen scapolottina: {filename}')
                        chosen_file_path = filename
                        control=1

    return chosen_file_path


def read_json(chosen_file_path:str):

    data = 0
    if os.path.isfile(chosen_file_path):
        with open(chosen_file_path, "r") as file:
            data = json.load(file)
    else:
        raise FileNotFoundError('File not found.')

    return data