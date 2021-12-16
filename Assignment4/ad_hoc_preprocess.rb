=begin

 - Author:
Alvaro Martinez Petit

 - Description:
This is an extra script for the fourth assignment of the Bioinformatics Programming Challenges Course. It performs some
transformations to the provided data in order for my mainScript.rb to function, despite being completely independent.
Execute the script by doing: $ ruby mainScript.rb arab_gene_file proteome_file

 @arab_gene_file = File name or file path of the arabidopsis CDS
 @proteome_file = File name or file path of the other proteome
 
After the execution, the following files will be created in the current directory:
 - TAIR10_proteome.fa
 - s_pombe_db.phr
 - s_pombe_db.pin
 - s_pombe_db.psq
 - tair10_db.phr
 - tair10_db.pin
 - tair10_db.psq

=end

require 'bio'

#########################################################################################################################
# WARNINGS
#########################################################################################################################

abort_msg = "The file 'TAIR10_proteome.fa' is going to be overwritten. "\
            "Please, move or remove this file from the current directory and try again."
abort(abort_msg) if File.file?('TAIR10_proteome.fa')

#########################################################################################################################
# INPUT CHECKING
#########################################################################################################################

arab_cds_file, proteome_file = ARGV

unless arab_cds_file && proteome_file
  abort "Please, run this script using the command 'ruby mainScript.rb arab_gene_file proteome_file' \n\n"
end

if ARGV.length > 2
    puts "\nINFO: this script only takes as input the two first positional arguments: 'arab_gene_file'and 'proteome_file', "\
          "in that specific order. Any other argument provided was ignored.\n\n"
end

abort_msg = "Please, provide a valid file path for the 'arab_gene_file' argument and try again. '#{arab_cds_file}' is not valid"
abort(abort_msg) unless File.file?(arab_cds_file)

abort_msg = "Please, provide a valid file path for the 'proteome_file' argument and try again. '#{proteome_file}' is not valid"
abort(abort_msg) unless File.file?(proteome_file)

#########################################################################################################################
# MAIN CODE
#########################################################################################################################

######## TRANSALATING TAIR10 CDS TO AMINOACIDS ##########

puts "Reading and translating arabidposis file..."

arab_data = Bio::FlatFile.auto(arab_cds_file)
arab_protein_data = File.new('TAIR10_proteome.fa', 'w')

arab_data = Bio::FlatFile.auto(arab_cds_file)
arab_protein_data = File.new('TAIR10_proteome.fa', 'w')

arab_data.each_entry do |entry|
  
  nucleotide_seq = entry.to_biosequence
  nucleotide_seq.na
  amino_seq = nucleotide_seq.translate
  
  arab_protein_data.write("\n>" + entry.definition )
  arab_protein_data.write("\n" + amino_seq)
  
end

arab_protein_data.close

puts "Translation completed"

######## CREATING BLAST DATABASES ##########

puts "Starting database creation..."

system "makeblastdb -in #{proteome_file} -dbtype 'prot' -out s_pombe_db"

system "makeblastdb -in TAIR10_proteome.fa -dbtype 'prot' -out tair10_db"

puts "Finished database creation"
