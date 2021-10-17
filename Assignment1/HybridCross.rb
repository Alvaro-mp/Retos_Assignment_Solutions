class HybridCross
  
  @@instances = []
  attr_accessor :parent_1
  attr_accessor :parent_2
  attr_accessor :f2_wild
  attr_accessor :f2_p1
  attr_accessor :f2_p2
  attr_accessor :f2_p1p2
  
  def initialize(params={})
    @parent_1 = params.fetch(:parent_1)
    @parent_2 = params.fetch(:parent_2)
    @f2_wild = params.fetch(:f2_wild, 0).to_i
    @f2_p1 = params.fetch(:f2_p1, 0).to_i
    @f2_p2 = params.fetch(:f2_p2, 0).to_i
    @f2_p1p2 = params.fetch(:f2_p1p2, 0).to_i
    
    @@instances << self
  end
  
  
  def HybridCross.instances
    return @@instances
  end
  
  
  def chi_square_test
    total_pop = f2_wild + f2_p1 + f2_p2 + f2_p1p2
    # expected populations:
    exp_wild = total_pop * (9.0/16.0)
    exp_p1 = total_pop * (3.0/16.0)
    exp_p2 = total_pop * (3.0/16.0)
    exp_p1p2 = total_pop * (1.0/16.0)
    
    result = ((f2_wild-exp_wild)**2)/exp_wild + ((f2_p1-exp_p1)**2)/exp_p1 + ((f2_p2-exp_p2)**2)/exp_p2 + ((f2_p1p2-exp_p1p2)**2)/exp_p1p2
    
    # as this problem has 3 degrees of freedom, significant p-values are found with a chi-square value above 7.815
    if result >= 7.815
      puts "Recording: #{parent_1.mutant_gene.name} is linked to #{parent_2.mutant_gene.name} with chi-square score of #{result}"
      parent_1.mutant_gene.linked_to << parent_2.mutant_gene.name
      parent_2.mutant_gene.linked_to << parent_1.mutant_gene.name
    end
  end
  
end