#!/usr/bin/env bash
#SBATCH --mem='40gb'
#SBATCH --cpus-per-task=5
#SBATCH --constraint=cal
#SBATCH --time='10:00:00'

hostname
source ~soft_bio_267/initializes/init_autoflow
source ~soft_bio_267/initializes/init_python
source ~soft_bio_267/initializes/init_degenes_hunter
source ~soft_bio_267/initializes/init_R

mode=$1
cur_dir=`pwd`
CODE_PATH=$cur_dir/scripts
datasets='aldh18a1;sqle;sprout;ischaemia'


if [ "$mode" == "0" ]; then
	mkdir final_counts
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/scRNA_SQLE/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/sqle_final_counts.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/scRNA_ALDH18A1/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/aldh18a1_final_counts.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/HUVEC/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/sprout_final_counts.txt
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/angiogenesis/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/ischaemia_final_counts.txt

	mkdir targets
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/scRNA_SQLE/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/sqle_target.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/scRNA_ALDH18A1/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/aldh18a1_target.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/HUVEC/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/sprout_target.txt
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/angiogenesis/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/ischaemia_target.txt
fi

if [ "$mode" == "1" ]; then
        echo "DEA, CEA and functional analysis"
        mkdir expResults
        packages='WDE'
        var_info=`echo -e "\\$packages=$packages,
		\\$final_counts_folder=$cur_dir/final_counts,
		\\$targets_folder=$cur_dir/targets,
                \\$datasets=$datasets" |  tr -d '[:space:]' `
        AutoFlow -w templates/deg_ceg_template.af -t '7-00:00:00' -m '40gb' -c 1 -o expResults -V $var_info -n 'cal' $2 #-u 1
fi

if [ "$mode" == "2" ]; then
        echo "Integrative analysis"
        mkdir intResults
        mkdir tmpResults
        mkdir hunterFolders
        mkdir plots
        dst=`echo $datasets | tr ";" " "`
        datasetsArray=( $dst )
     	# Para crear los simbolicos apoyarte en index_execution
        # for i in ${datasetsArray[@]}
        # do
        #         var=`grep $i $cur_dir/expResults/index_execution | cut -f 2`
        #         ln -s $var $cur_dir/hunterFolders/$i
        # done
        # Combinar resultados de la ejecución anterior.
        # De cada dataset buscar genes DEG en común y llevarlos a la plantilla
        $CODE_PATH/venndiagrams.R -n "`echo $cur_dir/hunterFolders/*/prevalent_degs/* | tr ' ' ","`" -o $cur_dir/plots/vennDiagram 
        get_venn_groups.rb -i "$cur_dir/hunterFolders/*/prevalent_degs/*" -o $cur_dir/tmpResults/group_assignation
        poss_genes=`cut -f 1 $cur_dir/tmpResults/group_assignation`
        grep ".*,.*,.*,.*" $cur_dir/tmpResults/group_assignation | cut -f 1 > $cur_dir/tmpResults/common_prev_genes
        # TODO PREPARAR GRÁFICAS VENN
	common_prev_genes=$cur_dir/tmpResults/common_prev_genes
        var_info=`echo -e "\\$hunterFolders=$cur_dir/hunterFolders,
                \\$datasets=$datasets,
                \\$poss_genes=$poss_genes,
                \\$common_prev_genes=$common_prev_genes" |  tr -d '[:space:]' `
        echo $var_info
        AutoFlow -w templates/integrativeAnalysis.af -t '7-00:00:00' -m '40gb' -c 1 -o 'intResults/'$datasets -V $var_info $2 -n 'cal' #-u 1
fi


if [ "$mode" == "3" ]; then
  	mkdir selected_PPIs
        mkdir datasets
        mkdir confident_results
        mkdir final_stats
        poss_genes=`cut -f 1 $cur_dir/tmpResults/group_assignation`
	prev_genes=$cur_dir/tmpResults/common_prev_genes
        # echo "Download and decompress STRING data (human)"
        # wget 'https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz' -O $cur_dir/datasets/string_data.txt.gz
        # gunzip $cur_dir/datasets/string_data.txt.gz
        # cut -d ' ' -f 1,2,16 $cur_dir/datasets/string_data.txt | tail -n +2 | sed 's/9606.//g' | tr ' ' "\t" > $cur_dir/datasets/human_ppt_interactions.txt
        # rm $cur_dir/datasets/string_data.txt

        echo 'Getting ENSP-gene Name relationships'

        # wget https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz -O $cur_dir/datasets/9606.protein.info.v11.5.txt.gz
        # gunzip $cur_dir/datasets/9606.protein.info.v11.5.txt.gz
        # cut -f 1,2 datasets/9606.protein.info.v11.5.txt | tail -n +2 > $cur_dir/datasets/ensp_geneName_dict.txt
        # sed -i 's/9606.//g' $cur_dir/datasets/ensp_geneName_dict.txt
        # rm $cur_dir/datasets/9606.protein.info.v11.5.txt
        
        echo 'Getting ENSG-ENSP relationships'

        # wget https://ftp.ensembl.org/pub/release-108/tsv/homo_sapiens/Homo_sapiens.GRCh38.108.uniprot.tsv.gz -O $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz
        # gunzip $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz
        # cut -f 1,3 $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv | tail -n +2 | awk '{print $2 "\t" $1}' > $cur_dir/datasets/ensp_ensg_dict.tsv
        # rm $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz

	echo "Prepare confident networks (PPI combined score >= 950) and select interactions"
        #awk '{if ($3 >= 950) print $0}' $cur_dir/datasets/human_ppt_interactions.txt > $cur_dir/tmpResults/confident_human_interactions.txt
        # Get all expressed genes
        #cut -f 1 $cur_dir/hunterFolders/*/*_DEA/filtered_count_data.txt | sort -u | tail -n +2 > $cur_dir/tmpResults/all_datasets_expressed_genes
        # Translate ENSG (expressed genes) to ENSP (to find matches in STRING interactions)
        # standard_name_replacer.py -I $cur_dir/datasets/ensp_ensg_dict.tsv -i $cur_dir/tmpResults/all_datasets_expressed_genes -c 1 -f 2 -t 1 -o $cur_dir/tmpResults/all_datasets_expressed_genes_transl -s "\t" 
        # Get blacklist (no-expressed genes)
        #grep -v -w -F -f $cur_dir/tmpResults/all_datasets_expressed_genes $cur_dir/hunterFolders/aldh18a1/aldh18a1_DEA/Common_results/* | cut -f 1 | tail -n +2 > $cur_dir/tmpResults/blacklist
        #standard_name_replacer.py -I $cur_dir/datasets/ensp_ensg_dict.tsv -i $cur_dir/tmpResults/blacklist -c 1 -f 2 -t 1 -o $cur_dir/tmpResults/blacklist_transl -s "\t"
        #grep -v -w -F -f $cur_dir/tmpResults/blacklist_transl $cur_dir/tmpResults/confident_human_interactions.txt > $cur_dir/tmpResults/confident_expressed_human_interactions.txt

        echo "Clustering networks (rber_pots) and selection"
        confident_ppi=$cur_dir/tmpResults/confident_expressed_human_interactions.txt
        cdlib_clusterize.py -i $confident_ppi -o confident_results/rber_pots.txt -m rber_pots -s -S confident_results/rber_pots_stats.txt
        #0       ENSG00000173692
        $CODE_PATH/get_metrics.rb -s string -i confident_results/rber_pots.txt -t clusters > final_stats/rber_pots_cluster_stats.txt
        # CON CLUSTERES HACER ENRIQUECIMIENTO FUNCIONAL
fi

exit


        exit
        

        #$CODE_PATH/code.rb -i
        

        #$CODE_PATH/string_ppi_selection.R -i $cur_dir/datasets/human_ppt_interactions.txt -c1 -c2 -k 950 -o $cur_dir/selected_PPIs/human_confident_PPI.txt

        # Select pairs in which both nodes are in the all datasaet expressed genes file
        #join by ENSP
        ###merge_tabular.py $cur_dir/datasets/ensp_ensg_dict.tsv $cur_dir/datasets/ensp_geneName_dict.txt > $cur_dir/tmpResults/ensp_ensg_geneName_dict 
        exit
        echo 'Translating STRING ENSP codes to ENSG'

        
        ###sed -i 's/\"//g' $cur_dir/datasets/human_ppi_ensembl_genes.txt
         echo 'Translating all expressed gene ENSG identifiers to gene name identifiers'
        # Translate all_datasets_expressed_genes ENSEMBL gene ID to HGNC gene ID:
        #standard_name_replacer.py -i $cur_dir/tmpResults/all_datasets_expressed_genes -I $cur_dir/datasets/protein_dictionary.txt -c 1 -f 1 -t 3 -o $cur_dir/tmpResults/all_datasets_expressed_genes_transl
        exit


        ###
