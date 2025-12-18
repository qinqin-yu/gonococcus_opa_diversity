#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20M
#SBATCH -p shared
#SBATCH -t 0-00:10
#SBATCH -o entrez-direct.out
#SBATCH -e entrez-direct.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --account=grad_lab

# Define function
function func1 {
  esearch -db sra -query "$1" \
    | elink -target biosample \
    | efetch -format docsum \
    | xtract -pattern DocumentSummary -element Accession
}

# Export function
export -f func1

#Call function
cat grad_lab_complete_genomes_short_read_sra_ids.txt \
  | xargs -n1 -I{} bash -c '
    sra="$1"
    biosample=$(func1 "$sra")
    printf "%s\t%s\n" "$sra" "$biosample"
  ' _ {} \
  >sra_to_biosample.tsv
  
# cat grad_lab_complete_genomes_short_read_sra_ids.txt | xargs -I {} bash -c 'func1 "{}"'>sra_to_biosample.tsv