=begin

Author:
Alvaro Martinez Petit

Description:
This is the main script for the first assignment of the Bioinformatics Programming Challenges Course.
Execute the script by doing: $ ruby mainScript.rb
In order for it to work, the following files must be in the same folder as the mainScript.rb file:
 - Gene.rb
 - SeedStock.rb
 - SeedStockDatabase.rb
 - HybridCross.rb
 - gene_information.tsv
 - seed_stock_data.tsv
 - cross_data.tsv
After the execution, a file named 'new_stock_file.tsv' will be created in the same folder as the script

=end

require './Gene.rb'
require './SeedStock.rb'
require './SeedStockDatabase.rb'
require './HybridCross.rb'

# ----------- TASK 1 ----------------------

# Loading seed stock data
stock_db = SeedStockDatabase.new
stock_db.load_from_files('./seed_stock_data.tsv', './gene_information.tsv')

# Planting 7 grams of each seed stock
stock_db.instances.each {|seed_stock|
  seed_stock.plant_grams(7)
}

# Saving the new database state to a new file
stock_db.write_database('./new_stock_file.tsv')

# ----------- TASK 2 --------------------

# Loading cross data
cross_file = File.new('./cross_data.tsv', 'r')
cross_data = cross_file.read.split("\n")
  # removing header
cross_data.delete_at(0)
cross_data.each {|line|
  HybridCross.new(
    parent_1: stock_db.get_seed_stock(line.split("\t")[0]),
    parent_2: stock_db.get_seed_stock(line.split("\t")[1]),
    f2_wild: line.split("\t")[2],
    f2_p1: line.split("\t")[3],
    f2_p2: line.split("\t")[4],
    f2_p1p2: line.split("\t")[5]
  )
}

# Perform chi-square test over all HybridCross instances
HybridCross.instances.each {|cross|
  cross.chi_square_test
}

puts
puts
puts "Final report:"
puts

stock_db.instances.each {|seed_stock|
  unless seed_stock.mutant_gene.linked_to.empty?
    puts "#{seed_stock.mutant_gene.name} is linked to #{seed_stock.mutant_gene.linked_to.join(", ")}"
  end
  
}

