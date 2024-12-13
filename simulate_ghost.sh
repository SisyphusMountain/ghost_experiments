#!/bin/bash
# by Theo Tricou

python Zombi/Zombi.py T SpeciesTreeParameters_no_ghost.tsv sim_no_ghost
python Zombi/Zombi.py T SpeciesTreeParameters_clade.tsv sim_clade
python Zombi/Zombi.py T SpeciesTreeParameters_sample.tsv sim_sample

python Zombi/Zombi.py G GenomeParameters.tsv sim_no_ghost
python Zombi/Zombi.py G GenomeParameters.tsv sim_clade
python Zombi/Zombi.py G GenomeParameters.tsv sim_sample

python Zombi/SpeciesSampler.py n 30 sim_no_ghost
python Zombi/SpeciesSampler.py i to_keep_clade.txt sim_clade
python Zombi/SpeciesSampler.py n 30 sim_sample

for i in sim_*; do
  cd "$i/SAMPLE_1"

  # Apply sed command to each .nwk file
  for file in *.nwk; do
    [ -f "$file" ] || continue  # Skip if no matching files
    sed -E -i \
      -e "s/([nR][a-zA-Z0-9]+)_[^:;]+/\\1/g" \
      -e "s/Root/9999/g" \
      -e "s/n([0-9]+)/100\\1/g" \
      "$file"
  done
  cd ../T
  if [ -f "CompleteTree.nwk" ]; then
    sed -E -i \
      -e "s/([nR][a-zA-Z0-9]+)_[^:;]+/\\1/g" \
      -e "s/Root/9999/g" \
      -e "s/n([0-9]+)/100\\1/g" \
      "CompleteTree.nwk"
  fi
  cd ../SAMPLE_1
  ALEobserve SampledSpeciesTree.nwk
  ALEml_undated SampledSpeciesTree.nwk SampledSpeciesTree.nwk.ale output_species_tree=y
  cp SampledSpeciesTree.nwk_SampledSpeciesTree.nwk.ale.spTree aletree
  ls *_sampledtree.nwk | parallel -P 0 ALEobserve
  ls *_sampledtree.nwk.ale | parallel -P 0 ALEml_undated aletree {}
  ls *uTs | parallel -P 0 'sed -i -e "s/([0-9]*[0-9])//g" -e "s/^\t//g" -e "/^#/d"'
  python3 ../../nodes_mappings_rec.py stat_ale_ghost_$i 
  cd ../..
done

Rscript ./plot_ghost.R

# GNU Ghost
