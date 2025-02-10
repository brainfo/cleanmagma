import gzip
import os
import pandas as pd
import shutil
import re
import logging

root_folder = '/mnt/data/hong/2022/DHJ1_human_obesity_placenta/data/gwas/egg-consortium.org'
## set logging path, and logging level
logging.basicConfig(filename=f'{root_folder}/mk_db.log', level=logging.DEBUG)
target_folder = '/mnt/storage/hong/2024/egg-consortium'


def unzip_file(file_path, target_folder):
    """unzip the gz file and move the gz file to a target folder"""
    file_name = file_path.replace('.gz', '')
    with gzip.open(file_path,"rb") as f_in, open(file_name,"wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
        ## move the zip file to a target folder
        shutil.move(file_path, os.path.join(target_folder, file_path.split("/")[-1]))

def mk_pval_file(file_path, n=42212):
    """
    return a file with fields rs_id, p, and n.
    """
    ## read the file, delimiter is tab
    df = pd.read_csv(file_path, sep='\t')
    ## get the column names
    col_names = df.columns.tolist()
    rsid_cols = [col for col in col_names if col.lower() in ['hm_rsid', 'rs_id','rsid','snp', 'markername', 'variant_id']]
    p_cols = [col for col in col_names if col.lower() in ['p', 'p-value', 'pval', 'p.value', 'p_value']]
    n_cols = [col for col in col_names if col.lower() in ['n', 'nsamples', 'totalsamplesize']]
    ## check if the rsid_cols, p_cols, and n_cols are unique
    if len(rsid_cols) != 1 or len(p_cols) != 1:
        ## not raise error, but print a warning and continue
        logging.warning(f"Warning: There should be only one rsid_col, p_col, and n_col in {file_path}")
        ## -[] what if not unique?
        return
    ## get the column names
    rsid_col = rsid_cols[0]
    p_col = p_cols[0]
    ## if n_cols is empty, use 42212
    ## simplified ifelse statement
    n_col = n_cols[0] if len(n_cols) > 0 else n
    ## create df2use with the column names ['rsid', 'p', 'n']
    df2use = pd.DataFrame({"rsid": df[rsid_col], "p": df[p_col], "n": df[n_col]})
    df2use['p'] = pd.to_numeric(df2use['p'], errors='coerce')
    df2use['n'] = pd.to_numeric(df2use['n'], errors='coerce')
    ## remove any nan row
    df2use = df2use[~df2use.isna().any(axis=1)]
    ## specify the column types using pandas for dataframe but not per column
    df2use = df2use.astype({'rsid': 'str', 'p': 'float', 'n': 'int'})
    ## write the file
    filename = f'{file_path.split("/")[-1]}_p.txt'
    df2use.to_csv(filename, index=False, sep='\t')

## for os walk 
# -[] maybe modulize to handle the breaks of the chain of programs

def mk_loc_file(file_path):
    """
    return a file with fields rs_id, chromosome,position.
    """
        ## read the file, delimiter is tab
    df = pd.read_csv(file_path, sep='\t')
    ## get the column names
    col_names = df.columns.tolist()
    rsid_cols = [col for col in col_names if col.lower() in ['rs_id','rsid','snp', 'markername', 'variant_id']]
    chromosome_cols = [col for col in col_names if col.lower() in ['chromosome', 'chrom', 'chromosome_name', 'chr']]
    ## using any(substring in string for substring in substring_list)
    position_cols = [col for col in col_names if any(substring in col.lower() for substring in ['position', 'pos', 'bp', 'basepair', 'base_pair_location', 'base_pair_location_start'])]
    if len(rsid_cols) != 1 or len(chromosome_cols) != 1 or len(position_cols) != 1:
        ## not raise error, but print a warning and continue
        logging.warning(f"Warning: There should be only one rsid_col, chromosome_col, and position_col in {file_path}")
        ## -[] what if not unique?
        return
    ## get the column names
    rsid_col = rsid_cols[0]
    chromosome_col = chromosome_cols[0]
    position_col = position_cols[0]
    ## create df2use with the column names ['rsid', 'chromosome', 'position']
    df2use = pd.DataFrame({"rsid": df[rsid_col], "chromosome": df[chromosome_col], "position": df[position_col]})
    ## coerce conversion for each column, pd.to_numeric, then astype
    df2use['chromosome'] = pd.to_numeric(df2use['chromosome'], errors='coerce')
    df2use['position'] = pd.to_numeric(df2use['position'], errors='coerce')
    ## remove any nan row
    df2use = df2use[~df2use.isna().any(axis=1)]
    df2use = df2use.astype({'rsid': 'str', 'chromosome': 'int', 'position': 'int'})
    filename = f'{file_path.split("/")[-1]}_loc.txt'
    df2use.to_csv(filename, index=False, sep='\t')

def run_single_gwas(root, file, target_folder, n=42212):
    unzip_file(os.path.join(root, file), target_folder)
    mk_pval_file(os.path.join(root, file.replace('.gz', '')), n)
    mk_loc_file(os.path.join(root, file.replace('.gz', '')))
    ## remove the unzipped file
    os.remove(os.path.join(root, file.replace('.gz', '')))

if __name__ == "__main__":
    for root, dirs, files in os.walk(root_folder):
        for file in files:
            if file.endswith('.gz'):
                ## logging
                logging.info(f"Processing {os.path.join(root, file)}")
                run_single_gwas(root, file)


