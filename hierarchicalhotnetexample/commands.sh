#!/usr/bin/env bash
git clone https://github.com/raphael-group/hierarchical-hotnet.git

data=$PWD/data
intermediate=$PWD/intermediate
results=$PWD/results

num_permutations=1000

################################################################################
#
#   Prepare data.
#
################################################################################

# Create intermediate data/results and results directories.
mkdir -p $intermediate
mkdir -p $results

for network in reactomefi2016 
do
    mkdir -p $intermediate/"$network"
done

for network in reactomefi2016
do
    for score in triglyceride_european 
    do
        mkdir -p $intermediate/"$network"_"$score"
    done
done

python scripts/format_score.py $data/Triglyceride.gene.stats.txt $data/triglyceride_european.tsv

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

echo "Construct similarity matrices..."

for network in reactomefi2016
do
    python src/construct_similarity_matrix.py \
        -i   $data/"$network"_edge_list.tsv \
        -o   $intermediate/"$network"/similarity_matrix.h5 \
        -bof $intermediate/"$network"/beta.txt
done

################################################################################
#
#   Permute data.
#
################################################################################

echo "Permuting scores..."

for network in reactomefi2016
do
    for score in triglyceride_european 
    do
        cp $data/"$score".tsv $intermediate/"$network"_"$score"/scores_0.tsv

        python src/find_permutation_bins.py \
            -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
            -igf $data/"$network"_index_gene.tsv \
            -elf $data/"$network"_edge_list.tsv \
            -ms  1000 \
            -o   $intermediate/"$network"_"$score"/score_bins.tsv

        for i in `seq $num_permutations`
        do
            python src/permute_scores.py \
                -i  $intermediate/"$network"_"$score"/scores_0.tsv \
                -bf $intermediate/"$network"_"$score"/score_bins.tsv \
                -s  "$i" \
                -o  $intermediate/"$network"_"$score"/scores_"$i".tsv
        done
    done
done

################################################################################
#
#   Construct hierarchies.
#
################################################################################

echo "Constructing hierarchies..."

for network in reactomefi2016 
do
    for score in triglyceride_european
    do
        for i in `seq 0 $num_permutations`
        do
            python src/construct_hierarchy.py \
                -smf  $intermediate/"$network"/similarity_matrix.h5 \
                -igf  $data/"$network"_index_gene.tsv \
                -gsf  $intermediate/"$network"_"$score"/scores_"$i".tsv \
                -helf $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv \
                -higf $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv
        done
    done
done

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

for network in reactomefi2016
do
    for score in triglyceride_european 
    do
        python src/process_hierarchies.py \
            -oelf $intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv \
            -oigf $intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_edge_list_"$i".tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/"$network"_"$score"/hierarchy_index_gene_"$i".tsv "; done) \
            -lsb  1 \
            -cf   $results/clusters_"$network"_"$score".tsv \
            -pl   $network $score \
            -pf   $results/sizes_"$network"_"$score".pdf
    done
done
