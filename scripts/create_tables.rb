#! /usr/bin/env ruby
#
## Dale el path entre comillas
## Script to create tables for HTML report
require 'optparse'
#
#
## METHODS

def load_module_genes(file)
	clusters_genes = {}
	File.open(file).each do |line|
		line.chomp!
		cluster, gene = line.split("\t")
		query = clusters_genes[cluster]
		if query.nil?
			clusters_genes[cluster] = [gene]
		else
			query << gene
		end
	end
	return clusters_genes
end


def save_table(output_path, table)
	File.open(output_path, 'w') do |f|
		table.each do |id, data|
			f.puts "#{id}\t#{data.join("\t")}"
		end
	end
end

def save_cluster_table(output_path, table)
	File.open(output_path, 'w') do |f|
		table.each do |id, data|
			data.each do |gene_info|
				f.puts "#{id}\t#{gene_info.join("\t")}"
			end
		end
	end
end

def save_intersecting_table(table, output_path)
	File.open(output_path, 'w') do |f|
		table.each do |data|
			f.puts "#{data.join("\t")}"
		end
	end
end

############################################################################################
## OPTPARSE
############################################################################################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  options[:cluster_files] = nil
  opts.on("-c", "--cluster_files PATH", "Path to STRING clustering files") do |item|
    options[:cluster_files] = item.split(",")
  end

  options[:intersecting_degs] = 'intersecting_degs.txt'
  opts.on("-d", "--intersecting_degs PATH", "Output path for intersecting STRING clusters genes and DEGs") do |item|
    options[:intersecting_degs] = item
  end

  options[:intersecting_cegs] = 'intersecting_cegs.txt'
  opts.on("-e", "--intersecting_cegs PATH", "Output path for intersecting STRING clusters genes and CEGs from selected modules") do |item|
    options[:intersecting_cegs] = item
  end

  options[:study_genes] = nil
  opts.on("-g", "--study_genes PATH", "Path to file with genes for each study") do |item|
    options[:study_genes] = item
  end

  options[:input_files] = nil
  opts.on("-i", "--input_files PATH", "Path to files for getting the different metrics") do |item|
    options[:input_files] = item.split(",")
  end

  options[:modules_list] = []
  opts.on("-l", "--modules_list STRING", "List of CEGs clusters selected. Example: sqle=35;sprout=13") do |item|
    options[:modules_list] = item.split(';').map{|a|a.split('=')}
  end

  options[:clustering_method] = 'rber'
  opts.on("-m", "--clustering_method STRING", "Clustering method to perform intersection analysis. Please choose: cpm, louvain, leiden, rber, rb") do |item|
    options[:clustering_method] = item
  end

  options[:network_type] = 950
  opts.on("-n", "--network_type STRING", "STRING network type to perform intersection analysis. Please choose: 950 (confident) or 300 (full)") do |item|
    options[:network_type] = item
  end

  options[:output_file] = 'output.txt'
  opts.on("-o", "--output_file PATH", "Output path") do |item|
    options[:output_file] = item
  end

  options[:output_matches] = 'output_matches.txt'
  opts.on("-O", "--output_matches PATH", "Output path with genes list matches in STRING clustering files") do |item|
    options[:output_matches] = item
  end

  options[:clustering_path] = nil
  opts.on("-p", "--clustering_path PATH", "Path to STRING clustering file") do |item|
    options[:clustering_path] = item
  end

end.parse!  

############################################################################################
## MAIN
#############################################################################################MAIN

experiment = {}
experiment_degs = {} # experiment => [genes]
experiment_cegs = {} # [experiment] = {clusterID => [genes]}

Dir.glob(File.join(options[:input_files], 'hunter_results_table_annotated.txt')).each do |file|
	experiment_name = file.split('/')[-3]
	expressed_genes = 0
	possible_deg = 0
	prevalent_deg = 0
	cegs_by_cluster = {}
	File.open(file).each do |line|
		line.chomp!
		next if line.include?('input_IDs')
		expressed_genes += 1
		fields = line.split("\t")
		ensembl_geneID = fields.first
		clusterID = fields[17]
		gene_tag = fields[15]
		if gene_tag == 'PREVALENT_DEG' || gene_tag == 'POSSIBLE_DEG'
			possible_deg += 1
		end
		if gene_tag == 'PREVALENT_DEG' 
			prevalent_deg += 1 
			query = experiment_degs[experiment_name]
			if query.nil?
				experiment_degs[experiment_name] = [ensembl_geneID]
			else
				query << ensembl_geneID
			end
		end
		query = cegs_by_cluster[clusterID]
		if query.nil?
			cegs_by_cluster[clusterID] = [ensembl_geneID]
		else
			query << ensembl_geneID
		end
	end
	experiment[experiment_name] = [expressed_genes - 1, possible_deg, prevalent_deg]
	experiment_cegs[experiment_name] = cegs_by_cluster
end

options[:modules_list].each do |experiment, moduleID|
	cegs = experiment_cegs[experiment]
	cegs.select!{|k,v| moduleID == k}
end


save_table(options[:output_file], experiment)

study_genes = []
File.open(options[:study_genes]).each do |line|
	study_genes << line.chomp!
end

string_clust_net_data = {}
degs_overlap = {}
cegs_overlap = {}
network_genes = []
experiment_genes_clusters = {}
cegs_by_cluster = {}


Dir.glob(File.join(options[:cluster_files], 'cluster_gene.txt')).each do |file|
	network_type = file.split('/')[-3].split('_').first
	clustering_method = file.split('/')[-3].split('_')[1]
	clusters_genes = load_module_genes(file)
	saved_data = []
	clusters_genes.each do |cluster, genes|
		common_genes = study_genes & genes
		cluster_genes = [cluster, genes]
		if common_genes.length > 0
			common_genes.each do |common_gene|
				saved_data << [common_gene, cluster, genes.length]
				if network_type == options[:network_type] && clustering_method == options[:clustering_method]
					experiment_genes_clusters[common_gene] = cluster_genes
				end
			end
		end
	end
	string_clust_net_data[network_type] = saved_data
end

save_cluster_table(options[:output_matches], string_clust_net_data)

intersecting_degs_table = []
intersecting_cegs_table = []

experiment_genes_clusters.each do |study_gene, cluster_data|
	cluster_ID, cluster_genes = cluster_data
	overlaps = []
	experiment_degs.each do |experiment, degs|
		intersecting_genes = (cluster_genes & degs).length
		overlaps << (intersecting_genes.fdiv(cluster_genes.length) * 100).round(3)
	end
	m_overlaps = []
	experiment_cegs.each do |experiment, module_data|
		module_data.each do |m_id, m_genes|
			intersecting_genes = (cluster_genes & m_genes).length
			m_overlaps << (intersecting_genes.fdiv(cluster_genes.length) * 100).round(3)
		end
	end
	intersecting_degs_table << [cluster_ID, cluster_genes.length, study_gene].concat(overlaps)
	intersecting_cegs_table << [cluster_ID, cluster_genes.length, study_gene].concat(m_overlaps)
end

save_intersecting_table(intersecting_degs_table, options[:intersecting_degs])
save_intersecting_table(intersecting_cegs_table, options[:intersecting_cegs])