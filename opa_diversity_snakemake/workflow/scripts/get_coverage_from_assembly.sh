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
  esearch -db assembly -query "$1" \
    | esummary \
    | xtract -pattern DocumentSummary -element Coverage
}

# Export function
export -f func1

#Call function
cat public_complete_genomes_assembly_accession.txt \
  | xargs -n1 -I{} bash -c '
    accession="$1"
    coverage=$(func1 "$accession")
    printf "%s\t%s\n" "$accession" "$coverage"
  ' _ {} \
  >public_complete_genomes_assembly_accession_to_coverage.tsv