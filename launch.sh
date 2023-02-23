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
datasets='aldh18a1;ischaemia;sprout;sqle'
dst=`echo $datasets | tr ";" " "`
datasetsArray=( $dst )
clustering_methods=( 'rber_pots' 'louvain' 'cpm' 'rb_pots' 'leiden' )

mkdir final_counts
mkdir targets
mkdir expResults
mkdir selected_PPIs
mkdir datasets
mkdir clustering_results
mkdir intResults
mkdir tmpResults
mkdir hunterFolders
mkdir plots

if [ "$mode" == "0" ]; then
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/scRNA_SQLE/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/sqle_final_counts.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/scRNA_ALDH18A1/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/aldh18a1_final_counts.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/HUVEC/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/sprout_final_counts.txt
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/angiogenesis/analysis/results/DEGenesHunter_results/ctrl_vs_mut/final_counts.txt $cur_dir/final_counts/ischaemia_final_counts.txt

	cp /mnt/home/users/bio_267_uma/beammzgz/projects/scRNA_SQLE/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/sqle_target.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/scRNA_ALDH18A1/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/aldh18a1_target.txt
	cp /mnt/scratch/users/bio_267_uma/beammzgz/HUVEC/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/sprout_target.txt
	cp /mnt/home/users/bio_267_uma/beammzgz/projects/angiogenesis/analysis/DEG_workflow_ori/TARGETS/ctrl_vs_mut_target.txt $cur_dir/targets/ischaemia_target.txt
fi

if [ "$mode" == "1" ]; then
        echo "DEA, CEA and functional analysis"
        packages='WDE'
        var_info=`echo -e "\\$packages=$packages,
		\\$final_counts_folder=$cur_dir/final_counts,
		\\$targets_folder=$cur_dir/targets,
                \\$datasets=$datasets" |  tr -d '[:space:]' `
        AutoFlow -w templates/deg_ceg_template.af -t '7-00:00:00' -m '40gb' -c 1 -o expResults -V $var_info -n 'cal' $2 #-u 1
fi

if [ "$mode" == "2" ]; then
        echo "Integrative analysis"

        ## Symbolic links to execution path:
        # for i in ${datasetsArray[@]}
        # do
        #         var=`grep $i $cur_dir/expResults/index_execution | cut -f 2`
        #         ln -s $var $cur_dir/hunterFolders/$i
        # done

        $CODE_PATH/venndiagrams.R -c -n "`echo $cur_dir/hunterFolders/*/prevalent_degs/* | tr ' ' ","`" -t "`ls hunterFolders/*/prevalent_degs/* | cut -d ' ' -f 9 | cut -d '/' -f 2 | tr "\n" ','`" -o $cur_dir/plots/vennDiagram_DEG 
        get_venn_groups.rb -i "$cur_dir/hunterFolders/*/prevalent_degs/*" -o $cur_dir/tmpResults/group_assignation
        cut -f 1 $cur_dir/tmpResults/group_assignation > $cur_dir/tmpResults/poss_genes
        grep ".*,.*,.*,.*" $cur_dir/tmpResults/group_assignation | cut -f 1 > $cur_dir/tmpResults/common_prev_genes
        
        # TODO: prepare variable with genes by study (the ones they characterised in Rohlenova, Rosano and Harvey)

        var_info=`echo -e "\\$hunterFolders=$cur_dir/hunterFolders,
                \\$datasets=$datasets,
                \\$common_prev_genes=$cur_dir/tmpResults/common_prev_genes" |  tr -d '[:space:]' `
        AutoFlow -w templates/integrativeAnalysis.af -t '7-00:00:00' -m '40gb' -c 1 -o intResults -V $var_info $2 -n 'cal' #-u 1       
fi


if [ "$mode" == "D" ]; then
        echo "Downloading datasets"
        poss_genes=`cut -f 1 $cur_dir/tmpResults/group_assignation`
	prev_genes=$cur_dir/tmpResults/common_prev_genes
        echo "Download and decompress STRING data (human)"
        wget 'https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz' -O $cur_dir/datasets/string_data.txt.gz
        gunzip $cur_dir/datasets/string_data.txt.gz
        cut -d ' ' -f 1,2,16 $cur_dir/datasets/string_data.txt | tail -n +2 | sed 's/9606.//g' | tr ' ' "\t" > $cur_dir/datasets/human_ppt_interactions.txt
        rm $cur_dir/datasets/string_data.txt

        echo 'Getting ENSP-gene Name relationships'

        wget https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz -O $cur_dir/datasets/9606.protein.info.v11.5.txt.gz
        gunzip $cur_dir/datasets/9606.protein.info.v11.5.txt.gz
        cut -f 1,2 datasets/9606.protein.info.v11.5.txt | tail -n +2 > $cur_dir/datasets/ensp_geneName_dict.txt
        sed -i 's/9606.//g' $cur_dir/datasets/ensp_geneName_dict.txt
        rm $cur_dir/datasets/9606.protein.info.v11.5.txt
        
        echo 'Getting ENSG-ENSP relationships'

        #wget https://ftp.ensembl.org/pub/release-108/tsv/homo_sapiens/Homo_sapiens.GRCh38.108.uniprot.tsv.gz -O $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz
        #gunzip $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz
        wget "https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz" -O $cur_dir/datasets/9606.protein.aliases.v11.5.txt.gz 
        gunzip $cur_dir/datasets/9606.protein.aliases.v11.5.txt.gz
        grep 'Ensembl_gene' $cur_dir/datasets/9606.protein.aliases.v11.5.txt | cut -f 1,2 | cut -d '.' -f 2 > $cur_dir/datasets/ensp_ensg_dict.tsv
        # cut -f 1,3 $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv | tail -n +2 | awk '{print $2 "\t" $1}' > $cur_dir/datasets/ensp_ensg_dict.tsv
        # rm $cur_dir/datasets/Homo_sapiens.GRCh38.108.uniprot.tsv.gz

        echo "Get all expressed genes"
        cut -f 1 $cur_dir/hunterFolders/*/*_DEA/filtered_count_data.txt | sort -u | tail -n +2 > $cur_dir/tmpResults/all_datasets_expressed_genes
        # Translate ENSG (expressed genes) to ENSP (to find matches in STRING interactions)
        standard_name_replacer.py -I $cur_dir/datasets/ensp_ensg_dict.tsv -i $cur_dir/tmpResults/all_datasets_expressed_genes -c 1 -f 2 -t 1 -o $cur_dir/tmpResults/all_datasets_expressed_genes_transl -s "\t" 
        # Get blacklist (no-expressed genes)
        grep -v -w -F -f $cur_dir/tmpResults/all_datasets_expressed_genes $cur_dir/hunterFolders/aldh18a1/aldh18a1_DEA/Common_results/* | cut -f 1 | tail -n +2 > $cur_dir/tmpResults/blacklist
        standard_name_replacer.py -I $cur_dir/datasets/ensp_ensg_dict.tsv -i $cur_dir/tmpResults/blacklist -c 1 -f 2 -t 1 -o $cur_dir/tmpResults/blacklist_transl -s "\t"
fi

if [ "$mode" == "3" ]; then
        echo "Prepare networks and select interactions"
        combScores=( 300 950 )
        for m in ${clustering_methods[@]}
        do
                for i in ${combScores[@]}
                do
                var_info=`echo -e "\\$network=$cur_dir/datasets/human_ppt_interactions.txt,
                        \\$common_prev_genes=$cur_dir/tmpResults/common_prev_genes,
                        \\$combScores=$i,
                        \\$clustering_methods=$m,
                        \\$dictionary=$cur_dir/datasets/ensp_ensg_dict.tsv,
                        \\$CODE_PATH=$CODE_PATH,
                        \\$cores=16,
                        \\$blacklist=$cur_dir/tmpResults/blacklist_transl" |  tr -d '[:space:]' `
                AutoFlow -w templates/clusteringAnalysis.af -t '7-00:00:00' -m '40gb' -c 1 -o clustering_results/"$i"_"$m" -V $var_info $2 -n 'cal' #-u 1
                done
        done
fi

if [ "$mode" == "4" ]; then
        source ~soft_bio_267/initializes/init_report_html
        echo "Generate result tables and reports"
        rm -rf final_stats
        mkdir final_stats
        path=`ls clustering_results/*/cdlib_clusterize.py_0000/cluster_stats.txt | cut -d ' ' -f 9 | tr "\n" ' '`
        merge_tabular.py $path > final_stats/all_methods_stats.txt
        sed -i "1i Param\t300_cpm\t300_leiden\t300_louvain\t300_rber_pots\t300_rb_pots\t950_cpm\t950_leiden\t950_louvain\t950_rber_pots\t950_rb_pots" final_stats/all_methods_stats.txt
        cluster_genes_path="$cur_dir/clustering_results/*_rber_pots/cdlib_clusterize.py_0000"
        $CODE_PATH/create_tables.rb -i "$cur_dir/hunterFolders/*/*_FEnr" -o final_stats/genes_expression_stats.txt -g $cur_dir/study_genes -c "$cluster_genes_path" -O $cur_dir/final_stats/output_matches.txt -n 950 -m rber -d $cur_dir/final_stats/degs_table.txt -e $cur_dir/final_stats/cegs_table.txt -l 'sqle=7;aldh18a1=34;sprout=30;ischaemia=30'
        $CODE_PATH/venndiagrams.R -c -n "`echo $cur_dir/intResults/cut_000*/coexp_genes/* | tr ' ' ","`" -t "`ls $cur_dir/intResults/cut_000*/coexp_genes/* | cut -d '/' -f 14| tr "\n" ','`" -o $cur_dir/plots/vennDiagram_CEG
        $CODE_PATH/venndiagrams.R -c -n "`echo $cur_dir/tmpResults/*_CEG_DEG | tr ' ' ","`" -t "`ls $cur_dir/intResults/cut_000*/coexp_genes/* | cut -d '/' -f 14| tr "\n" ','`" -o $cur_dir/plots/vennDiagram_DEG_CEG
        report_html -t templates/html_template.erb -d final_stats/genes_expression_stats.txt,final_stats/all_methods_stats.txt,final_stats/output_matches.txt,final_stats/cegs_table.txt,final_stats/degs_table.txt
fi
