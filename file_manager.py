import os
import shutil
import logging

def gather_gff_files(source_dirs, target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    for src_dir in source_dirs:
        for file in os.listdir(src_dir):
            if file.endswith('.gff'):
                src_file_path = os.path.join(src_dir, file)
                target_file_path = os.path.join(target_dir, file)
                if not os.path.exists(target_file_path):
                    shutil.move(src_file_path, target_file_path)
                logging.info(f"Přemístěno {src_file_path} do {target_file_path}")
