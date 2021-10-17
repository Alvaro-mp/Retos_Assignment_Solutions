class SeedStockDatabase
  
  attr_accessor :instances
  
  def initialize
    @instances = []
  end  
  
  
  def load_from_files(stockfilepath, genefilepath)
    stockfile = File.new(stockfilepath, "r")
    genefile = File.new(genefilepath, "r")
    
    genes = genefile.read
    stocks = stockfile.read.split("\n")
    # removing header
    stocks.delete_at(0)
    
    stocks.each {|seed_stock|
      # if the mutant gene referenced in the Seed Stock file is present in the Gene file
      coincidence = Regexp.new(/(#{seed_stock.split("\t")[1]})\t(.+)\t(.+)/).match(genes)
      if coincidence 
        gene = Gene.new(
          id: coincidence[1],
          name: coincidence[2],
          phenotype: coincidence[3]
        )
        
        @instances.push(SeedStock.new(
          id: seed_stock.split("\t")[0],
          mutant_gene: gene,
          last_planted: seed_stock.split("\t")[2],
          storage: seed_stock.split("\t")[3],
          grams_remaining: seed_stock.split("\t")[4]
        ))
      end
        
    }
    
    genefile.close
    stockfile.close
  end
  
  
  def get_seed_stock(id)
    instances.each {|stock|
      if stock.id == id
        return stock
      end
      }
    puts "Seed stock with id \"#{id}\" not found in database"
    
  end
  
  
  def write_database(filepath)
    file = File.new(filepath, "w")
    file.write("Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining")
    instances.each {|stock|
      file.write("\n#{stock.id}\t#{stock.mutant_gene.id}\t#{stock.last_planted}\t#{stock.storage}\t#{stock.grams_remaining}")
      
      }
    
    file.close
  end
  
end