=begin

Author:
Alvaro Martinez Petit

Description:
This is the main script for the third assignment of the Bioinformatics Programming Challenges Course.
Execute the script by doing: $ ruby mainScript.rb target_repeat gene_file
 @target_repeat = DNA sequence of the targeted repeat of the insertional mutagenesis. It must follow the IUPAC coding for DNA sequences
 @gene_file = File with a list of genes in which the 'target_repeat' will be searched. They must be all separated by the newline '\n' character
After the execution, the following files will be created (or overwritten) in the current directory:
 - question_4a_mutagenesis_target_local.gff
 - genes_without_target.txt
 - question_5_mutagenesis_target_global.gff

=end

require 'bio'
require 'stringio'
require 'rest-client'

#########################################################################################################################
# WARNINGS
#########################################################################################################################

abort_msg = "The file 'question_4a_mutagenesis_target_local.gff' is going to be overwritten. "\
            "Please, move or remove this file from the current directory and try again."
abort(abort_msg) if File.file?('question_4a_mutagenesis_target_local.gff')

abort_msg = "The file 'genes_without_target.txt' is going to be overwritten. "\
            "Please, move or remove this file from the current directory and try again."
abort(abort_msg) if File.file?('genes_without_target.txt')

abort_msg = "The file 'question_5_mutagenesis_target_global.gff' is going to be overwritten. "\
            "Please, move or remove this file from the current directory and try again."
abort(abort_msg) if File.file?('question_5_mutagenesis_target_global.gff')

#########################################################################################################################
# INPUT CHECKING
#########################################################################################################################

search_seq, genefile = ARGV

unless search_seq && genefile
  abort "Please, run this script using the command 'ruby mainScript.rb target_repeat gene_file' \n\n"
end

if ARGV.length > 2
    puts "\nINFO: this script only takes as input the two first positional arguments: 'target_repeat' and 'gene_file', "\
          "in that specific order. Any other argument provided was ignored.\n\n"
end

search_seq.split("").each{|letter|
  unless "ACGTRYSWKMBDHVN.-".include?(letter.upcase)
    abort "Unexpected character #{letter} in 'target_repeat'.\n"\
          "Please, provide a valid 'target_repeat' following the IUPAC coding for DNA sequences."
  end
}

abort_msg = "Please, provide a valid file path for the 'gene_file' argument and try again. '#{genefile}' is not valid"
abort(abort_msg) unless File.file?(genefile)

#########################################################################################################################
# FUNCTIONS
#########################################################################################################################


def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end


def get_target_indices(string, bio_re)
  # inspired from 'https://stackoverflow.com/questions/5241653/ruby-regex-match-and-get-positions-of'
  return string.enum_for(:scan, bio_re).map{ [Regexp.last_match.begin(0), Regexp.last_match.end(1)-1].map{|n| n.to_i} }
end

# After using this function I found the Bio::Location class from BioRuby, which is much more powerful
def get_location_info(position)
  strand = "+"
  strand = "-" if position.include?("complement")
  
  re = Regexp.new(/(\d+)..(\d+)/)
  positions = position.scan(re).map{|match| match.map{|index| index.to_i}}
  
  return strand, positions
end


#########################################################################################################################
# MAIN CODE
#########################################################################################################################

# Load the content of the gene file and standardize gene case
gene_list = IO.readlines(genefile, chomp: true).map{|gene| gene.capitalize}

# Create a Bio::Sequence from the provided target sequence and create regular expresions for '+' and '-' strands
repeat = Bio::Sequence::NA.new(search_seq)
pattern = Regexp.new(/(?=(#{repeat.to_re}))/)
pattern_complement = Regexp.new(/(?=(#{repeat.complement.to_re}))/)

# Empty array to store entries of all genes in the list
entries = []

# Empty hash to store chormosomic gene locations and avoid fetching the database again later
gene_locations = Hash.new

puts "Starting sequence search"

gene_list.each do |gene|
  # DEBUG PRINT
  puts "Starting search for gene #{gene_list.index(gene)+1}/#{gene_list.length}"
  
  response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&id=#{gene}&format=embl")
  
  unless response
    puts "It was not possible to retrieve information for gene #{gene}. Skipping..."
    next
  end
  
  # Create a Bio::FlatFile providing EMBL format, as it was the format fetched, and a stream with the content
  datafile = Bio::FlatFile.new(Bio::EMBL, StringIO.new(response.body))
  datafile.each_entry do |entry|
    next unless entry.accession
    # Retrieve strand and coordinates of the current gene from the text in the entry description.
    # Strand information is useless and will be overwritten
    strand, positions = get_location_info(entry.description)
    # Add an entry to the gene_location Hash with the current gene as key and a simple Array of Integers as values
    gene_locations[gene] =  positions.flatten # if there is a join() in the description this fails
    
    # Create an editable Bio::Sequence from the entry in which add new Bio::Features
    bioseq = entry.to_biosequence
    # I iterate over the entry instead of the biosequence to avoid the risk of modify anything unintentionally
    entry.features.each do |feature|
      # Iterate only over exon features related to the current gene
      next unless feature.feature == "exon"
      next unless feature.assoc["note"].downcase.include?(gene.downcase)
      # Get strand and coordinates for the current exon relative to the gene
      strand, positions = get_location_info(feature.position)
      # In case the exon is split in various segments we need to search in each of them
      positions.each do |position|
        start_pos, end_pos = position
        # Get the exon sequence from the whole gene sequence
        exon_seq = entry.seq[start_pos-1..end_pos-1]
        # Get indices of all targets in the exon sequence.
        # If the strand is '-' we look for the target's complement in the positive strand, to avoid doing tedious math with indices
        if strand == "+"
          targ_indices = get_target_indices(exon_seq, pattern)
        else
          targ_indices = get_target_indices(exon_seq, pattern_complement)
        end
        
        # For each target matched...
        targ_indices.each do |match|
          # ...calculate coordinates relative to the gene...
          targ_positions = (start_pos + match[0]).to_s + ".." + (start_pos + match[1]).to_s
          targ_positions = "complement(" + targ_positions + ")" if strand == "-"
          # ...and create a new Bio::Feature, with some qualifiers, including the parent gene for later use
          new_feature = Bio::Feature.new('mutagenesis_target', targ_positions)
          new_feature.append(Bio::Feature::Qualifier.new('repeat_motif', repeat.seq.to_s))
          new_feature.append(Bio::Feature::Qualifier.new('type', 'nucleotide_motif'))
          new_feature.append(Bio::Feature::Qualifier.new('strand', strand))
          new_feature.append(Bio::Feature::Qualifier.new('parent', gene))
          # There might be exons overlapping, so we don't need to add features that already exist
          # Nonetheless, the sentence below doesn't prevent from adding repeated features as it compares objects themselves
          # and not their content, so I will remove repeated features later
          bioseq.features |= [new_feature]
        end      
      end
    end
    # Add the modified Bio::Sequence to the Array of entries
    entries << bioseq
  end
end

#################### CREATION OF GFF FILE #1 #########################################################

gff_file = File.new("question_4a_mutagenesis_target_local.gff", "w")
gff_file.write("##gff-version 3\n")

entries.each do |entry|
  entry.features.each do |feature|
    # Iterate only over the added features
    next unless feature.feature == "mutagenesis_target"
    # Get coordinates for the current feature
    strand, positions = get_location_info(feature.position)
    
    # Retrieve the number of chromosome from the gene AGI locus code
    seqid = feature.assoc["parent"][2]
    source = "."
    type = feature.assoc["type"]
    # There is no "join()" in the feature position for sure, so we can flatten the Array of Arrays
    start = positions.flatten[0]
    end_ = positions.flatten[1]
    score = "."
    strand = feature.assoc["strand"]
    phase = "."
    attributes = "Parent=#{feature.assoc["parent"]}"
    
    # Join all the values with '\t' characteres and add the new line
    line_content = [seqid, source, type, start, end_, score, strand, phase, attributes].join("\t")
    gff_file.write("\n#{line_content}")
  end
end

gff_file.close

# Remove the repeated entries from the file, as they are easier to handle as Strings
gff_input = File.new("question_4a_mutagenesis_target_local.gff", "r")
gff = gff_input.read()
gff_input.close

# I will write this content in a new file AND I will use it in the rest of the code
gff = gff.split("\n").uniq.join("\n")

gff_file = File.new("question_4a_mutagenesis_target_local.gff", "w")
gff_file.write(gff)
gff_file.close

#################### REPORT ON WHICH GENES HAVE THE REPEAT ########################################

# Create a GFF3 Object from the content of the GFF file created
repeats = Bio::GFF::GFF3.new(gff)

# Empty Array to store genes that contain the repeat
genes_with_repeat = []

# Retrieve the Parent gene on each line of the file and add it to the Array
repeats.records.each do |line|
  next unless line.feature == "nucleotide_motif"
  genes_with_repeat << line.get_attribute('Parent')
end

# Remove duplicated genes and standardize case
genes_with_repeat = genes_with_repeat.uniq.map{|gene| gene.capitalize}

# Find the genes present in the original Array and not present in the new one
genes_without_repeat = gene_list - genes_with_repeat

# Write the asked report
report = File.new("genes_without_target.txt", "w")
report.write("Provided genes that do NOT contain the specified target sequence:")
genes_without_repeat.each do |gene|
  report.write("\n· #{gene}")
end
# If there is no genes without the repeat, write "None"
report.write("\n\nNone") if genes_without_repeat.empty?
  
report.close

#################### CREATION OF GFF FILE #2 #########################################################

gff_file = File.new("question_5_mutagenesis_target_global.gff", "w")
gff_file.write("##gff-version 3\n")

# Iterate again over the lines of the first GFF file
repeats.records.each do |line|
  next unless line.feature == "nucleotide_motif"
  
  # Keep all values except for the start and end coordinates 
  seqid = line.seqid # I could not find what value should I give the seqid in this file different from the previous
  source = line.source
  type = line.feature_type
  # Use the gene chromosome coordinates retrieved at the beginning to calculate the feature chromosome coordiantes
  start = gene_locations[line.get_attribute('Parent')][0] + (line.start - 1)
  end_ = gene_locations[line.get_attribute('Parent')][0] + (line.end - 1)
  score = line.score
  strand = line.strand
  phase = line.phase
  # Don't keep the Parent information as it produces errors when uploading to EnsEMBL
  attributes = nil
  
  # Replace the nil values with "." and join the values with "\t" character, then write a new line
  line_content = [seqid, source, type, start, end_, score, strand, phase, attributes].map{|x| x||"."}.join("\t")
  gff_file.write("\n#{line_content}")
end

gff_file.close

puts "\n\nExecution finished!\n"