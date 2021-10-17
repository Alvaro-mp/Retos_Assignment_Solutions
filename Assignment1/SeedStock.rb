class SeedStock

  attr_accessor :id
  attr_accessor :mutant_gene
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining
  
  
  def initialize(params = {})
    @id = params.fetch(:id, "None")
    @mutant_gene = params.fetch(:mutant_gene)
    @last_planted = params.fetch(:last_planted, "00/00/0000")
    @storage = params.fetch(:storage, "not provided")
    @grams_remaining = params.fetch(:grams_remaining, 0).to_i
    # it would be better to cast to float as it is possible to plant fractions of grams,
    # but I casted to integer to keep the similitude with the original file
  end
  
  
  def plant_grams(grams)
    
    unless grams_remaining == 0
      today = Time.new
      @last_planted = "#{today.day}/#{today.month}/#{today.year}"
    end
    
    unless grams_remaining - grams <= 0
      @grams_remaining -= grams
    else
      @grams_remaining = 0
      puts "WARNING: we have run out of Seed Stock \"#{id}\""
    end
  end


end