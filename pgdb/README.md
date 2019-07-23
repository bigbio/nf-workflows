Protein Databases Generation using *pypgatk*
============================================

Run the workflow::
-------
	
	git clone https://github.com/bigbio/nf-workflows.git

	nextflow run main.nf --tool_basepath path_to_pygatk_base_dir --cosmic_user_name username --cosmic_password password 
 
 
Required tools:
---------

pypgatk: https://pgatk.readthedocs.io/en/latest/pypgatk.html#installation

gsutil: https://cloud.google.com/storage/docs/gsutil_install

git-lfs: https://github.com/git-lfs/git-lfs/wiki/Installation
