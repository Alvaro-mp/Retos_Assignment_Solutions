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

class InteractionNetwork
  
  @@instances = []
  attr_accessor :interactors
  attr_accessor :depth
  attr_accessor :cutoff_score
  attr_accessor :annotations
  
  
  def initialize(params={})
    @interactors = params.fetch(:interactors)
    @depth = params.fetch(:depth, nil)
    @cutoff_score = params.fetch(:cutoff_score, nil)
    @annotations = params.fetch(:anotations, [])
    
    @@instances << self
  end
  
  
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

  def annotate()
  
    go_response = fetch("www.geneontology.org/ontology/subsets/goslim_plant.obo")
    
    unless go_response
      puts "It was not possible to retrieve the Gene Ontology. Stop execution and retry"
    else
      gene_ontology = go_response.body
    end
    
    
    self.interactors.each{|gene|

      annotation = Annotation.new(AGI_locus_code: gene)

      uniprot_response = fetch("http://togows.org/entry/ebi-uniprot/#{gene}")

      if uniprot_response
        go_ids = uniprot_response.body.scan(/GO:[\d]{7}/)
        go_ids = go_ids.uniq
        /KEGG;(?<kegg_id>.*);/ =~ uniprot_response.body
        kegg_id = kegg_id.strip

        kegg_response = fetch("http://togows.org/entry/kegg-genes/#{kegg_id}.json")
        if kegg_response
          kegg_data = JSON.parse(kegg_response.body)
          pathways = kegg_data[0]["pathways"]
          pathways.keys.each{|key|
            annotation.add_KEGG(key, pathways[key])
            }

        else
          puts "It was not possible to retrieve KEGG information for gene #{gene}"
        end

        go_ids.each{|go_id|
          pattern = Regexp.new(/#{go_id}\nname:(?<term_name>.*)\nnamespace:(?<namespace>.*)\n/)
          if pattern.match?(gene_ontology)
            captures = pattern.match(gene_ontology).captures
            term_name = captures[0].strip
            namespace = captures[1].strip

            if namespace == "biological_process"
              annotation.add_GO(go_id, term_name)
            end
          end
          }

      else
        puts "It was not possible to retrieve uniprot information for gene #{gene}"

      end
      
      self.annotations << annotation
      }   

  end
  
  def print_network(file)
    
    genes = self.interactors.join(", ")
    kegg_annot = []
    go_annot = []

    self.annotations.each{|annotation|
      annotation.KEGG.each{|kegg| kegg_annot << kegg}
      annotation.GO.each{|go| go_annot << go}
      }
    
    kegg_annot = kegg_annot.uniq
    go_annot = go_annot.uniq
    
    file.write("Interactor genes: #{genes}\n")
    
    file.write("KEGG pathways involved:\n")
    kegg_annot.each{|k_annot|
      file.write("id: #{k_annot["id"]}, name:#{k_annot["pathway"]}\n")
      }
    
    file.write("Related GO terms:\n")
    go_annot.each{|g_annot|
      file.write("id: #{g_annot["id"]}, term:#{g_annot["term_name"]}\n")
      }
        
  end

  
  
end