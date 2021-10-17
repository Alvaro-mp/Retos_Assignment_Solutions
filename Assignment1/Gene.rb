class Gene
  
  attr_accessor :id
  attr_accessor :name
  attr_accessor :phenotype
  attr_accessor :linked_to
  
  
  def initialize(params = {})
    unchecked_id = params.fetch(:id, "None")
    gene_pattern = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/)
    
    # if the id is "None" or doesn't match the pattern, the code stops.
    unless gene_pattern.match(unchecked_id)
        ###raise Exception.new "Gene id \"#{unchecked_id}\" is not valid. Execution stopped."
        abort("Error: Gene id \"#{unchecked_id}\" is not valid. Execution stopped.")
    end
    
    @id = unchecked_id
    @name = params.fetch(:name, "not provided")
    @phenotype = params.fetch(:phenotype, "not provided")
    @linked_to = []
  end

end