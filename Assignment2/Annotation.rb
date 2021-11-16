# == Annotation
#
# This is a class to represent annotations for genes,
# including KEGG and GO identifiers
#
# == Summary
# 
# This is a class to create annotations and modify
# existing ones
#

class Annotation
  
  # Get/Set the identifier of the Arabidopsis gene
  # @!attribute [rw]
  # @return [String] AGI Locus Code
  attr_accessor :AGI_locus_code
  
  # Get/Set KEGG annotations
  # @!attribute [rw]
  # @return [Array] KEGG annotations as Hashes
  attr_accessor :KEGG
  
  # Get/Set GO annotations
  # @!attribute [rw]
  # @return [Array] GO annotations as Hashes
  attr_accessor :GO
  
  
  # Create a new instance of Annotations. Only the AGI is required

  # @param params [Hash] hash containing an AGI Locus Code as String, an Array of KEGG annotations
  # and an Array of GO annotations
  # @return [Annotation] an instance of Annotation
  def initialize(params={})
    @AGI_locus_code = params.fetch(:AGI_locus_code)
    @KEGG = params.fetch(:KEGG, [])
    @GO = params.fetch(:GO, [])
  end
  
  
   # Add a KEGG annotation to an existing Annotation object

  # @param id [String] KEGG identifier of the pathway referenced
  # @param pathway [String] Name of the KEGG pathway
  # @return 
  def add_KEGG(id, pathway)
    self.KEGG << {"id" => id, "pathway" => pathway}
  end
  
  
  # Add a GO annotation to an existing Annotation object

  # @param id [String] GO identifier of the term referenced
  # @param term_name [String] Name of the GO term
  # @return
  def add_GO(id, term)
    self.GO << {"id" => id, "term_name" => term}
  end
  
end