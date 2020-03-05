#!/bin/bash
#SBATCH --account=ctb-yxia
#SBATCH --time=1-00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/project/ctb-yxia/ghadie84/edgotype_fitness_effect_full_model/code/%x-%j.out
#SBATCH --job-name=ghadie84_produce_template_annotated_interactome_1

source ~/venv/bin/activate
python produce_template_annotated_interactome_1.py