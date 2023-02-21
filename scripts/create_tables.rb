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

  options[:study_genes] = nil
  opts.on("-g", "--study_genes PATH", "Path to file with genes for each study") do |item|
    options[:study_genes] = item
  end

  options[:input_files] = nil
  opts.on("-i", "--input_files PATH", "Path to files for getting the different metrics") do |item|
    options[:input_files] = item.split(",")
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

Dir.glob(File.join(options[:input_files], 'hunter_results_table_annotated.txt')).each do |file|
	experiment_name = file.split('/')[-3]
	expressed_genes = 0
	possible_deg = 0
	prevalent_deg = 0
	File.open(file).each do |line|
		line.chomp!
		expressed_genes += 1
		fields = line.split("\t")
		gene_tag = fields[15]
		possible_deg += 1 if gene_tag == 'PREVALENT_DEG' || gene_tag == 'POSSIBLE_DEG'
		prevalent_deg += 1 if gene_tag == 'PREVALENT_DEG'
	end
	experiment[experiment_name] = [expressed_genes - 1, possible_deg, prevalent_deg]
end

save_table(options[:output_file], experiment)

study_genes = []
File.open(options[:study_genes]).each do |line|
	study_genes << line.chomp!
end

string_clust_net_data = {}

Dir.glob(File.join(options[:cluster_files], 'cluster_gene.txt')).each do |file|
	network_type = file.split('/')[-3].split('_').first
	clusters_genes = load_module_genes(file)
	saved_data = []
	clusters_genes.each do |cluster, genes|
		common_genes = study_genes & genes
		if common_genes.length > 0
			common_genes.each do |common_gene|
				saved_data << [common_gene, cluster, genes.length]
			end
		end
	end
	string_clust_net_data[network_type] = saved_data
end


save_cluster_table(options[:output_matches], string_clust_net_data)
