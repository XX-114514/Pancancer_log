#!/bin/bash -l
#SBATCH --job-name=c2l-env
#SBATCH --partition=regular
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out
# Builds ${PC_ROOT}/envs/conda/c2l: cell2location 0.1.5 (scvi-tools>=1.3), squidpy, celltypist 1.7.1
# (same version as the atlas CellTypist labels). Separate from the frozen replay environments.
set -euo pipefail
set -a; source ~/.config/pancancer/pancancer.local.sh; set +a
module load micromamba
export MAMBA_ROOT_PREFIX=/vast/scratch/users/$USER/micromamba
P=$PC_ROOT/envs/conda/c2l
[ -e "$P" ] && { echo "exists: $P" >&2; exit 2; }
micromamba create -y -p "$P" -c conda-forge --override-channels python=3.11 pip
eval "$(micromamba shell hook -s bash)"; micromamba activate "$P"
pip install --no-cache-dir cell2location==0.1.5 squidpy celltypist==1.7.1 scikit-learn
# cell2location 0.1.5 imports pyro deep_getattr, removed after pyro-ppl 1.9.0
pip install --no-cache-dir pyro-ppl==1.9.0
pip freeze > "$P/pip_freeze.txt"
python -c "import cell2location,scvi,torch,squidpy,celltypist;print(cell2location.__version__,scvi.__version__,torch.__version__,torch.version.cuda,squidpy.__version__,celltypist.__version__)"
