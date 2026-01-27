# Adding Partial Charge Files

The provided `conda_lock.yml` has been generated for a `linux-64` machine. Future generation of future partial charges should also be run on the same architecture.

In a directory containing the ligand.sdf (without charges) should have a directiort names `_generate_charges_network` with a script like `run.sh`:

```bash
#!/bin/bash

group_name="jacs_set"
sys_name="bace"
nagl="am1bcc/openff-gnn-am1bcc-1.0.0.pt"

# Define the current directory and charge methods
current_dir="$(pwd)"
charge_methods=("am1bcc_at" "am1bccelf10_oe" "nagl_off" "am1bcc_oe")
conda_prefix="$(conda info --base)"
python_version="$(python --version | awk '{print $2}' | cut -d. -f1,2)"
nagl_models=(None None "${conda_prefix}/envs/openfe-benchmarks/lib/python${python_version}/site-packages/openff/nagl_models/models/${nagl}" None)

# Iterate through all charge methods
for i in "${!charge_methods[@]}"; do
    charge_method="${charge_methods[$i]}"
    nagl_model="${nagl_models[$i]}"

    echo "Using charge method: $charge_method"

    # Prepare nagl_model argument
    nagl_arg=""
    if [ "$nagl_model" != "None" ]; then
        nagl_arg="--nagl-model $nagl_model"
    fi

    # Run the charge molecules script
    python ../../../../data_generation/charge_molecules.py \
        --input-path ../ligands.sdf \
        --output-dir ../ \
        --charge-method $charge_method \
        $nagl_arg \
        --n-cores $SLURM_CPUS_ON_NODE > log_${charge_method}.txt
done

python ../../../../data_generation/generate_industry_lomap_networks.py \
    --system-group "${group_name}" \
    --system-name "$sys_name" \
    --input-sdf ../ligands.sdf \
    --out-dir ../ > log_network.txt

```
