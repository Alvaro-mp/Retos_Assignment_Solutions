=begin

 - Author:
Alvaro Martinez Petit

 - Description:
This is the main script for the fourth assignment of the Bioinformatics Programming Challenges Course.
Execute the script by doing: $ ruby mainScript.rb proteome_file_1 proteome_file_2 database_proteome_1 database_proteome_2

 @proteome_file_1 = File name or file path of the first proteome to search in
 @proteome_file_2 = File name or file path of the second proteome to search in
 @database_proteome_1 = File name or file path of the BLAST database containing only the sequences of the first proteome provided
 @database_proteome_2 = File name or file path of the BLAST database containing only the sequences of the second proteome provided
 
After the execution, the following files will be created  in the current directory:
 - ortholog_candidates.txt


=end

require 'bio'
require 'stringio'
require 'rest-client'

#########################################################################################################################
# FUNCTIONS
#########################################################################################################################


# Print the estimated time remaining of an iterative task

# @param start_time [Integer] 
# @param current_loop [Integer] 
# @param total_loops [Integer]
# @return [nil] prints a verbose sentence to the standard output
def print_time_remaining(start_time, current_loop, total_loops)
  remaining_loops = total_loops - current_loop
  current_time = Time.now.to_i
  seconds_remaining = (current_time - start_time) * remaining_loops / current_loop 
  hours = seconds_remaining / 3600
  minutes = (seconds_remaining - 3600*hours ) / 60
  seconds = seconds_remaining - 3600*hours - 60*minutes

  if hours > 0
    puts "Estimated finishing time: #{hours}h, #{minutes} minutes and #{seconds} seconds"
  else
    puts "Estimated finishing time: #{minutes} minutes and #{seconds} seconds"
  end
  
end

#########################################################################################################################
# WARNINGS
#########################################################################################################################

abort_msg = "The file 'ortholog_candidates.txt' is going to be overwritten. "\
            "Please, move or remove this file from the current directory and try again."
abort(abort_msg) if File.file?('ortholog_candidates.txt')

#########################################################################################################################
# INPUT CHECKING
#########################################################################################################################

proteome_1_file, proteome_2_file, proteome_1_db, proteome_2_db = ARGV

unless proteome_1_file && proteome_2_file && proteome_1_db && proteome_2_db
  abort "Please, run this script using the command 'ruby mainScript.rb proteome_file_1 proteome_file_2 database_proteome_1 database_proteome_2' \n\n"
end

if ARGV.length > 4
    puts "\nINFO: this script only takes as input the four first positional arguments: 'proteome_file_1', 'proteome_file_2', "\
          "'database_proteome_1' and 'database_proteome_2', in that specific order. Any other argument provided was ignored.\n\n"
end

abort_msg = "Please, provide a valid file path for the 'proteome_file_1' argument and try again. '#{proteome_1_file}' is not valid"
abort(abort_msg) unless File.file?(proteome_1_file)

abort_msg = "Please, provide a valid file path for the 'proteome_file_2' argument and try again. '#{proteome_2_file}' is not valid"
abort(abort_msg) unless File.file?(proteome_2_file)

#########################################################################################################################
# MAIN CODE
#########################################################################################################################

proteome_1_factory = Bio::Blast.local('blastp', proteome_1_db, '-b 1 -K 100')
proteome_2_factory = Bio::Blast.local('blastp', proteome_2_db, '-b 1 -K 100')
# to decide the selected BLAST paramenters, I looked for all the possible parameters in:
# http://etutorials.org/Misc/blast/Part+V+BLAST+Reference/Chapter+13.+NCBI-BLAST+Reference/13.3+blastall+Parameters/

# Hash to store ortholog candidate pairs
ortholog_candidates = Hash.new()

# In order to speed up the search, I need to find the sortest proteome to iterate over its sequences
file_1_seqs = 0
file_2_seqs = 0
File.open(proteome_1_file, "r") {|file|  file_1_seqs = file.read.scan(/^>/).size}
File.open(proteome_2_file, "r") {|file|  file_2_seqs = file.read.scan(/^>/).size}

# I also store the total number of sequences to output the search progress properly
if file_1_seqs < file_2_seqs
  shortest_proteome = proteome_1_file
  total_seqs = file_1_seqs
else
  shortest_proteome = proteome_2_file
  total_seqs = file_2_seqs
end

# I create two hashes to link files to each other and with their respective factory
opposite_file_to = {proteome_1_file => proteome_2_file, proteome_2_file => proteome_1_file}
factory_for_proteome = {proteome_1_file => proteome_1_factory, proteome_2_file => proteome_2_factory}

# I retrieve the file for the largest proteome using the previously declared hash
largest_proteome = opposite_file_to[shortest_proteome]

# SEARCHING INFO PRINT
progress_flag = 5
puts "\nStarting ortholog candidate search. Producing feedback every #{progress_flag}% of progress."

# ITERATIVE CODE
shortest_proteome_proteins = Bio::FlatFile.auto(shortest_proteome)

# I declare some variables to output the search progress properly
current_protein_nb = 0
starting_timestamp = Time.now.to_i

# Loop over each entry of the shortest proteome
shortest_proteome_proteins.each_entry do |protein|

  # DEBUG PRINT
  current_protein_nb += 1
  percentage = current_protein_nb*100.0/total_seqs
  
  if percentage%progress_flag < 100.0/total_seqs
    puts "Evaluated #{percentage.round()}% of proteome." 
    print_time_remaining(starting_timestamp, current_protein_nb, total_seqs)
  end

  # Get current protein sequence
  sequence = protein.seq

  # Perform the query in the appropriate database
  report = factory_for_proteome[largest_proteome].query(">myseq\n#{sequence.to_s}")

  unless report.hits().empty?
    
    # Get best hit
    best_hit = report.hits()[0]

    # In case the e-value of the best hit is significant, perform the reciprocal query
    if best_hit.evalue() < 0.01
      # Remove hyphens from the target sequence to avoid annoying warnings
      ungapped_target_seq = best_hit.target_seq.tr("-","")
      # Perform reciprocal query
      report = factory_for_proteome[shortest_proteome].query(">myseq\n#{ungapped_target_seq}")

      unless report.hits().empty?
        
        # Get the reciprocal best hit
        reciprocal_best_hit = report.hits()[0]

        # In case the reciprocal best hit equals the protein in the current loop, add the protein and its best hit to the ortholog candidate hash
        if reciprocal_best_hit.definition == protein.definition
          ortholog_candidates[protein.entry_id] = best_hit.hit_id
        end

      else
        puts "Strange behavior with protein #{protein.entry_id}. "\
          "Its best hit didn't return any hits."
      end
    end
  end

end

##OUTPUT A REPORT WITH CANDIDATES
report_file = File.new("ortholog_candidates.txt", "w")
report_file.write("Ortholog candidates found:\n")

ortholog_candidates.each do |key, value|
    report_file.write("\n· #{key} <=> #{value}")
end

report_file.write("\nNone") if ortholog_candidates.empty?

report_file.close

=begin

The next step in an orthology search would be look for common annotations to both candidates
This could be accomplished the same way we did in Assignment 2, retrieving KEGG and GO annotations
for the involved genes.

=end
