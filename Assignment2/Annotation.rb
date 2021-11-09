class Annotation
  
  attr_accessor :AGI_locus_code
  attr_accessor :KEGG
  attr_accessor :GO
  
  def initialize(params={})
    @AGI_locus_code = params.fetch(:AGI_locus_code)
    @KEGG = params.fetch(:KEGG, [])
    @GO = params.fetch(:GO, [])
  end
  
  def add_KEGG(id, pathway)
    self.KEGG << {"id" => id, "pathway" => pathway}
  end
  
  def add_GO(id, term)
    self.GO << {"id" => id, "term_name" => term}
  end
  
end