#!/usr/bin/env python
import os
import h5py


def find_hdf5_files(root_path, filename):
    """ Generator to find specific files in directory tree """
    for dirpath, dirnames, files in os.walk(root_path):
        if filename in files:
            yield os.path.join(dirpath, filename)

def repack_hdf5(original_file, temp_file):
    """ Repacks HDF5 file to reduce size by removing deleted spaces """
    with h5py.File(original_file, 'r') as old_h5, h5py.File(temp_file, 'w') as new_h5:
        for item in old_h5.keys():
            old_h5.copy(item, new_h5)
    os.rename(temp_file, original_file)  # Replace original file with the repacked version
    print(f"Repacked and replaced original file at {original_file}")

# Path to the root directory where the search will start
root_path = '/opt/conda/'
file_name = 'Dfam.h5'

# Find the file and delete the 'Families' group
for file_path in find_hdf5_files(root_path, file_name):
    print(f"Found HDF5 file at {file_path}")
    with h5py.File(file_path, 'r+') as file:
        if 'Families' in file:
            del file['Families']
            print("'Families' group deleted from:", file_path)

            # Specify a temporary file path for the repacked file
            temp_file_path = file_path.replace('.h5', '_temp.h5')
            repack_hdf5(file_path, temp_file_path)
        else:
            print("'Families' group does not exist in the file at:", file_path)
            # repack anyway
            repack_hdf5(file_path, temp_file_path)




# try to delete unnecessary files:
to_delete = ['rpsblast', 'deltablast', 'rpstblastn',
             'psiblast', 'tblastn','blastb', 'tblastx',
             'blastx']


for fname in to_delete:
    for file_path in find_hdf5_files(root_path, fname):
        print(f"Found file {fname} at {file_path}")
        os.remove(file_path)

