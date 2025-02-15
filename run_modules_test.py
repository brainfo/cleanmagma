import mk_db

root_folder = '/mnt/data/hong/2022/DHJ1_human_obesity_placenta/data/gwas/egg-consortium.org'
## set logging path, and logging level
target_folder = '/mnt/storage/hong/2024/egg-consortium'

mk_db.run_single_gwas(root_folder, 'Fetal_Effect_European_meta_NG2019.txt.gz', target_folder, p = ['p'], n=['n_ownbw'], rsid=['rsid'], chromosome=['chr'], position=['pos'])