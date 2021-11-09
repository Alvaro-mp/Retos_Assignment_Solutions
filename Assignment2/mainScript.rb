require './InteractionNetwork.rb'
require 'yaml'

gene_list = IO.readlines("ArabidopsisSubNetwork_GeneList.txt", chomp: true)

interaction_networks = InteractionNetwork.search_networks(gene_list = gene_list,
                                                          max_depth = 3,
                                                          min_intact_score = 0.5)

interaction_networks.each {|network|
    network.annotate()
  }

File.open("report.txt", "w") {|out_file|
    interaction_networks.each {|network|
            out_file.write("New interaction network\n\n")
            network.print_network(out_file)
            out_file.write("\n\n")
    }
}
