#! /usr/bin/env ruby
require 'optparse'

############################################################################################
## METHODS
############################################################################################

def load_tabular_file(file, skip = 0)
  records = []
  File.open(file).each do |line|
    line.chomp!
    fields = line.split("\t")
    records << fields 
  end
  records.shift(skip) unless skip == 0
  return records
end

def load_value(hash_to_load, key, value, unique = true)
 	query = hash_to_load[key]
  if query.nil?
      value = [value] if value.class != Array
      hash_to_load[key] = value
  else
      if value.class == Array
          query.concat(value)
      else
          query << value
      end
      query.uniq! unless unique == nil
  end
end

def median(array)
  return nil if array.empty?
  sorted = array.sort
  len = sorted.length
  (sorted[(len - 1) / 2] + sorted[len / 2]) / 2.0
end


def build_metrics(filename, filedata, type, sample, source, expanded, assoc_path = nil, pairs_data = nil)
	pairs_num = filedata.length
	nodes = filedata.flatten.uniq
	phen_nodes = nodes.to_s.scan(/HP/).length
	net_basename = sample.split("_").first	
	if type == 'net'
		hp = nodes.count{|n| n.match(/^HP:/)}
		dis_nodes = nodes.length - hp
		if filename.include?("filtered_HP")
			puts ["#{sample}_net", source, expanded, "raw_d2p_pairs", pairs_num].join("\t")
			puts ["#{sample}_net", source, expanded, "raw_d2p_Pnodes", phen_nodes].join("\t")
			puts ["#{sample}_net", source, expanded, "raw_d2p_Dnodes", dis_nodes].join("\t")	
		else
			puts ["#{sample}_net", source, expanded, "cleaned_d2p_pairs", pairs_num].join("\t")
			puts ["#{sample}_net", source, expanded, "cleaned_d2p_Pnodes", phen_nodes].join("\t")
			puts ["#{sample}_net", source, expanded, "cleaned_d2p_Dnodes", dis_nodes].join("\t")	
		end
	elsif type == 'assoc'
		if filename.include?("hp")
			puts ["#{sample}_net", source, expanded, "raw_p2p_assoc", pairs_num].join("\t")
			puts ["#{sample}_net", source, expanded, "raw_p2p_nodes", phen_nodes].join("\t")
		elsif filename.include?("no_cluster")
			puts ["#{sample}_net", source, expanded, "no_cluster_pairs", pairs_num].join("\t")	
		elsif filename.include?("reliable_pairs") && !filename.include?("cleaned_reliable_pairs")
			puts ["#{sample}_net", source, expanded, "thr_p2p_assoc", pairs_num].join("\t")
			puts ["#{sample}_net", source, expanded, "thr_p2p_nodes", phen_nodes].join("\t")
		elsif filename.include?("cleaned_reliable_pairs")
			puts ["#{sample}_net", source, expanded, "nr_p2p_assoc", pairs_num].join("\t")
			puts ["#{sample}_net", source, expanded, "nr_p2p_nodes", phen_nodes].join("\t")	
		end	
	elsif type == 'clusters'
		clusters = {}
		filedata.each do |cluster, hp|
			load_value(clusters, cluster, hp)
		end
		cluster_sizes = clusters.values.map {|hps| hps.length}
		cluster_nodes = clusters.values.flatten
		unique_nodes = cluster_nodes.uniq
		nodes_count = cluster_nodes.group_by(&:itself).values.map{|x| x.length}
		mean_nodes = nodes_count.inject(0, :+).fdiv(nodes_count.length)
		mean_size = cluster_sizes.inject(0, :+).fdiv(cluster_sizes.length)

		all_nodes_hash = unique_nodes.group_by(&:itself)
   	clust_pairs = pairs_data.select{|nodeA, nodeB| all_nodes_hash.include?(nodeA) && all_nodes_hash.include?(nodeB)}
	
		puts ["cl_num", clusters.keys.length].join("\t")
		puts ["cl_size_min", cluster_sizes.min].join("\t")
		puts ["cl_size_max", cluster_sizes.max].join("\t")
		puts ["cl_size_median", median(cluster_sizes)].join("\t")
		puts ["cl_size_mean", mean_size].join("\t")
		puts ["cl_nodes_mean", mean_nodes].join("\t")	
		puts ["cl_nodes", unique_nodes.length].join("\t")
		puts ["cl_pairs", clust_pairs.length].join("\t")
	end		
end

############################################################################################
## OPTPARSE
############################################################################################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  options[:sample] = nil
  opts.on("-s", "--sample PATH", "Sample used to perform metrics") do |item|
    options[:sample] = item
  end

  options[:input_files] = nil
  opts.on("-i", "--input_files PATH", "Path to files for getting the different metrics") do |item|
    options[:input_files] = item.split(",")
  end

  options[:type] = nil
  opts.on("-t", "--type STRING", "Type of data to be analyzed") do |item|
    options[:type] = item
  end

  options[:assoc_path] = nil
  opts.on("-a", "--assoc_path PATH", "Path to the last filtered assoc file") do |item|
    options[:assoc_path] = item
  end

  options[:pairs_file] = nil
  opts.on("-p", "--pairs_file PATH", "Path to file containing the pairs between nodes, used for calculating total pairs in cluster metrics") do |item|
    options[:pairs_file] = item
  end

end.parse!  

############################################################################################
## MAIN
############################################################################################
pairs_data = []
pairs_data = load_tabular_file(options[:pairs_file]) if !options[:pairs_file].nil?
pairs_data.map!{|pair| [pair[0], pair[1]]}

source = []
expanded = nil
options[:input_files].each do |file|
	if !File.exist?(file)
		warn("WARN: File #{file} not exists")
		next
	end
	data = load_tabular_file(file)
	filename = File.basename(file)
	if options[:sample].include?("MONDO")
		sampl_source = options[:sample].split("2").first	
		source = "#{sampl_source}\tyes"
	else
		source = "#{options[:sample]}\tno"
	end
	if options[:sample].include?("exp")
		expanded = "yes"
	else
		expanded = "no"
	end		
	build_metrics(filename, data, options[:type], options[:sample], source, expanded, options[:method], pairs_data)
end	

if options[:assoc_path]
	puts ["#{options[:sample]}_net", source, expanded, "assoc_path", options[:assoc_path]].join("\t")
end		