require 'json'
require 'stringio'
require 'rest-client'
require './Annotation.rb'

# Create a function called "fetch" that we can re-use everywhere in our code

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


# == InteractionNetwork
#
# This is a class to represent InteractionNetworks through
# its interactors, depth of the search, cutoff interaction
# score and annotations
#
# == Summary
# 
# This is a class to generate gene Interaction Networks 
# and handle different aspects of them, as annotate or
# output them
#

class InteractionNetwork
  
  # Class variable containing all instances
  # @!attribute [rw]
  # @return [Array] All instances
  @@instances = []
  
  # Get/Set the interactors array
  # @!attribute [rw]
  # @return [Array] All interactors of the InteractionNetwork
  attr_accessor :interactors
  
  # Get/Set the network's search depth
  # @!attribute [rw]
  # @return [Integer] search depth
  attr_accessor :depth
  
  # Get/Set the network's minimum score for interactions
  # @!attribute [rw]
  # @return [Float] The minimum intact-score
  attr_accessor :cutoff_score
  
  # Get/Set network's annotations
  # @!attribute [rw]
  # @return [Array] The annotations for all interactors in the network
  attr_accessor :annotations
  
  
    
  # Create a new instance of InteractionNetwork

  # @param params [Hash] hash containing an Array of interactors as String, the search depth as Integer,
  # the minimum intact-score as Float and an Array of Annotations
  # @return [InteractionNetwork] an instance of InteractionNetwork
  def initialize(params={})
    @interactors = params.fetch(:interactors)
    @depth = params.fetch(:depth, nil)
    @cutoff_score = params.fetch(:cutoff_score, nil)
    @annotations = params.fetch(:anotations, [])
    
    @@instances << self
  end
  
  
  # Function that searches for networks in a list of provided AGI Locus Codes, according to the provided maximum depth
  # and minimum intact-score to consider the interaction valid. It has the limitation of not checking previously
  # initialized networks; the function only creates new ones from the list of genes

  # @param gene_list [Array] gene pool in which is needed to look for interactors, as String
  # @param max_depth [Integer] the maximum depth for the search
  # @param min_intact_score [Float] the minimum intact-score to consider an interaction as valid
  # @return [Array] array of all InteractionNetworks created
  def self.search_networks(gene_list, max_depth, min_intact_score)
    
    # There is no point in searching networks for a single gene, so in that case the program stops and warns the user
    unless gene_list.is_a?(Array)
      abort("Network Search Error: Please, provide a list of AGI Locus Codes to search links between them")
    end
    
    # to ensure consistency I format all provided genes
    gene_list = gene_list.map{|gene| gene.capitalize()}
    
    # Array to store interaction networks
    networks = []
    
    gene_list.each{|gene|
      
      # DEBUG PRINT
      puts "Starting search for gene #{gene_list.index(gene)+1}/#{gene_list.length}"
      
      # do a recursive interactor search for each gene in the list, specifiyng a maximum depth for the search and
      # a minimum score to consider the interaction as valid. In addition, an array with all the genes in upper
      # levels is passed to avoid doing a recursive search on genes that are going to have a deeper search. 
      interactor_list = recursive_search(gene, max_depth, min_intact_score, current_depth=1, upper_genes=gene_list)
      
      # the 'recursive_search' method provides capitalized strings as output, so there is no need to format it
      ### interactor_list = interactor_list.map{|interactor| interactor.capitalize()}
      
      networks.each{|network|
        
        # if current network and the interactor list have at least one element in common, merge the two arrays and
        # delete the old network from the array of networks
        unless (network & interactor_list).empty?
          interactor_list |= network
          networks.delete_at(networks.index(network)) if networks.include?(network)
        end 
        }
      
      networks << interactor_list
      
      }
    
    # I only want to find the interactor genes from the original list, so I transform each network to the intersection
    # of itself and the original list
    networks = networks.map{|network| network & gene_list}
    
    interaction_networks = []
    
    # finally, I create an InteractionNetwork object for each interaction network in the list. Those InteractionNetwork
    # objects are stored either into an array to be returned as into their class instances array
    networks.each{|network|
      
      inter_network = InteractionNetwork.new(
        interactors: network,
        depth: max_depth,
        cutoff_score: min_intact_score,
        )
      
      interaction_networks << inter_network
      }
    
    return interaction_networks
  
  end


  # Function to call on a recursive way to search for interactors of a provided gene until the maximum depth is reached

  # @param gene [String] the code of the gene as a String
  # @param max_depth [Integer] the maximum depth for the search as a Integer
  # @param min_intact_score [Float] the minimum intact-score to consider an interaction as valid, as a Float
  # @param current_depth [Integer] the current depth of the search as a Integer
  # @param upper_genes [Array] array of String containing the codes of genes in upper levels of the search
  # @param found_interactors [Array] array of String containing the codes of previously found genes that interact with the current
  # @return [Array] array of String containing the current and previously found interactors
  def self.recursive_search(gene, max_depth, min_intact_score, current_depth=0, upper_genes=[], found_interactors=[])
    
    # retrieve a list of interactor genes according to some databases
    gene_interactors = get_interactors(gene, min_intact_score)
  
    # Before merging arrays into each other, it is safe to ensure format consistency in their elements
    gene = gene.capitalize()
    upper_genes = upper_genes.map{|upper| upper.capitalize()}
    gene_interactors = gene_interactors.map{|interactor| interactor.capitalize()}
    
    # In case the array of interactors is empty, the own gene is added
    gene_interactors |= [gene]
  
    # in case the depth limit is not reached
    if current_depth < max_depth
        
      gene_interactors.each{|interactor|
        
        # if the current interactor is present in any of the upper levels, a recursive interactor search will go
        # deeper from that level, so there is no point in doing it at this point. Otherwise, do a recursive search
        unless upper_genes.include?(interactor)
          
          # in the next recursive search, the current_depth is increased in 1 unit and the current interactors are 
          # now upper_genes too
          found_interactors |= recursive_search(gene = interactor, 
                                                max_depth = max_depth,
                                                min_intact_score = min_intact_score, 
                                                current_depth = current_depth+1, 
                                                upper_genes = (upper_genes|gene_interactors), 
                                                found_interactors = found_interactors)
          
        end      
        }
      
    end
    
    # return the union of the previously found interactors and the currently found interactors, so that there is no
    # duplicated elements
    return (found_interactors | gene_interactors)
    
  end


  # funtion that retrieves the available interactors in a set of databases for a specific gene according to a minimum
  # intact-score

  # @param gene [String] the code of the gene as a String
  # @param min_intact_score [Float] the minimum intact-score to consider an interaction as valid, as a Float
  # @return [Array] array of String containing the found interactors for the current gene
  def self.get_interactors(gene, min_intact_score)
    
    interactors = []
    
    # the search size could be expanded at this point, including other REST services that output files in the
    # 'tab25' format
    databases = [
      "http://bar.utoronto.ca:9090/psicquic/webservices/current/search/interactor/#{gene}?format=tab25",
      "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{gene}?format=tab25"
      ]
    
    databases.each{|url|
  
      response = fetch(url)
  
      if response
  
        content = StringIO.new response.body
  
        # until we reach the end of file...
        until content.eof?
  
          # retrieve one line at a time
          line = content.gets
          /intact-miscore:(?<score>\d\.?\d*)/ =~ line
  
          unless score.to_f < min_intact_score
            # get all AGI Locus Codes metioned in the line
            interactors |= line.scan(/A[Tt]\d[Gg]\d\d\d\d\d/)
          end
  
        end
  
      else
        puts "Response error when trying to retrieve interactors for #{gene} from #{url}. Interactors not retrieved." 
      end
    }
    
    return interactors
    
  end


  # function that annotates the genes of a InteractionNetwork with the associated KEGG pathways and GO terms
  
  # @return
  def annotate()
  
    # retrieve the Gene Ontology
    go_response = fetch("www.geneontology.org/ontology/go.obo")
    
    unless go_response
      puts "It was not possible to retrieve the Gene Ontology. Stop execution and retry"
    else
      gene_ontology = go_response.body
    end
    
    
    self.interactors.each{|gene|

      # instanciate new Annotation for current gene
      annotation = Annotation.new(AGI_locus_code: gene)

      uniprot_response = fetch("http://togows.org/entry/ebi-uniprot/#{gene}")

      if uniprot_response
        # get all GO ids without duplicates
        go_ids = uniprot_response.body.scan(/GO:[\d]{7}/)
        go_ids = go_ids.uniq
        # get the KEGG id
        /KEGG;(?<kegg_id>.*);/ =~ uniprot_response.body
        kegg_id = kegg_id.strip

        # get all pathways from the KEGG database for the current gene
        kegg_response = fetch("http://togows.org/entry/kegg-genes/#{kegg_id}.json")
        if kegg_response
          kegg_data = JSON.parse(kegg_response.body)
          pathways = kegg_data[0]["pathways"]
          pathways.keys.each{|key|
            # add current pathway to the annotation of the current gene
            annotation.add_KEGG(key, pathways[key])
            }

        else
          puts "It was not possible to retrieve KEGG information for gene #{gene}"
        end

        go_ids.each{|go_id|
          # pattern to find GO ids in the Gene Ontology and capture the required data
          pattern = Regexp.new(/#{go_id}\nname:(?<term_name>.*)\nnamespace:(?<namespace>.*)\n/)
          if pattern.match?(gene_ontology)
            captures = pattern.match(gene_ontology).captures
            term_name = captures[0].strip
            namespace = captures[1].strip

            # only add GO term if it belongs to the 'biological process' part
            if namespace == "biological_process"
              annotation.add_GO(go_id, term_name)
            end
          end
          }
        
        # After doing all this code I found an easier way to retrieve just the GO:IDs and GO:term_names
        # from the ebi-uniprot response without requesting the whole Gene Ontology, a way that also included
        # requesting the ebi-uniprot data in JSON format, but I don't have the time to implement it and the
        # solution above works too. Nonetheless, my solution would be useful if we wanted to extract any other
        # information from GO entries apart from IDs, Term names and Namespace.

      else
        puts "It was not possible to retrieve uniprot information for gene #{gene}"

      end
      
      self.annotations << annotation
      }   

  end
  
  
  # function that outputs a report about a specific InteractionNetwork
  
  # @param file [File/String] Either an open File object or a file path as String
  # in which write the content of the InteractionNetwork
  # @return 
  def print_network(file)
    
    # open the File if a String is provided, setting a variable to close it
    to_close = false
    if file.is_a?(String)
      file = File.new(file, "w")
      to_close = true
    end
    
    if file.is_a?(File)
      
      # variables that store the values to print
      genes = self.interactors.join(", ")
      kegg_annot = []
      go_annot = []
  
      # store all the annotations of the same database in the same container
      self.annotations.each{|annotation|
        annotation.KEGG.each{|kegg| kegg_annot << kegg}
        annotation.GO.each{|go| go_annot << go}
        }
      
      # remove duplicates in annotations
      kegg_annot = kegg_annot.uniq
      go_annot = go_annot.uniq
      
      # print data
      file.write("Interactor genes: #{genes}\n")
      
      file.write("KEGG pathways involved:\n")
      kegg_annot.each{|k_annot|
        file.write(" - id: #{k_annot["id"]}, name: #{k_annot["pathway"]}\n")
        }
      
      file.write("Related GO terms:\n")
      go_annot.each{|g_annot|
        file.write(" + id: #{g_annot["id"]}, term: #{g_annot["term_name"]}\n")
        }
      
      # if the file was open inside this function, close it
      if to_close
        file.close
      end
    
    # in case the parameter provided is neither a String nor a File, warn the user
    else
      abort("Argument error: Please provide either the file path as String or an open File object.")
    end
    
  end

  
  
end