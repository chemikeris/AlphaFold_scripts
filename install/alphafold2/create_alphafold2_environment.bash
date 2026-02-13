alphafold_src_directory='/data/alphafold'
conda env create -f ./conda_alphafold2.yml
conda activate alphafold2
cd $alphafold_src_directory
pip install -r requirements.txt
pip3 install --upgrade --no-cache-dir jax==0.4.26 jaxlib==0.4.26+cuda12.cudnn89  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
python -c "import jax; print(jax.devices())"
echo 'to get working comit, use this command: git checkout cd0357af72396da3b75fe4c9cf86f8d9ba4e4d4e'
cd -
