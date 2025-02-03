import os
import shutil

def make_dirs():
    dirs_created = []
    # Get the current working directory
    current_dir = os.getcwd()
    
    # List all files in the current directory
    files_in_directory = [f for f in os.listdir(current_dir) if os.path.isfile(f)]
    
    # Filter files that contain 'pre-dictionary' in their name
    pre_dict_files = [f for f in files_in_directory if 'pre-dictionary' in f]
    
    # Loop through each file with 'pre-dictionary' in the name
    for file in pre_dict_files:
        # Grab the part after the first underscore
        parts = file.split('_', 1)  # Split into two parts at the first underscore
        if len(parts) > 1:
            subpart = parts[1].split('.')[0] + '.' # Get the part after the first underscore and before the extension
            
            # Create the directory if it doesn't already exist
            subpart_dir = os.path.join(current_dir, subpart)
            if not os.path.exists(subpart_dir):
                os.makedirs(subpart_dir)
                dirs_created.append(subpart)
            # Loop through all files again and move those containing the subpart into the directory
            for other_file in files_in_directory:
                if subpart in other_file:  
                    source = os.path.join(current_dir, other_file)
                    destination = os.path.join(subpart_dir, other_file)
                    shutil.move(source, destination)
    return dirs_created
