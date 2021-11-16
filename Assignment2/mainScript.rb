=begin

Author:
Alvaro Martinez Petit

Description:
This is the main script for the second assignment of the Bioinformatics Programming Challenges Course.
Execute the script by doing: $ ruby mainScript.rb
In order for it to work, the following files must be in the same folder as the mainScript.rb file:
 - Annotation.rb
 - InteractionNetwork.rb
 - ArabidopsisSubNetwork_GeneList.txt
After the execution, a file named 'report.txt' will be created in the same folder as the script

=end

require './InteractionNetwork.rb'

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
