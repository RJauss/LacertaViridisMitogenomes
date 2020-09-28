#!/usr/bin/env ruby
# frozen_string_literal: true

require 'bio'
require 'bio/db/nexus'
require 'fileutils'

base_dir = 'MITOS_2019_08_09'

gen_hash = {
  'nad1' => {
    size: nil
  },
  'nad2' => {
    size: nil
  },
  'cox1' => {
    size: nil
  },
  'cox2' => {
    size: nil
  },
  'atp8' => {
    size: nil
  },
  'atp6' => {
    size: nil
  },
  'cox3' => {
    size: nil
  },
  'nad3' => {
    size: nil
  },
  'nad4l' => {
    size: nil
  },
  'nad4' => {
    size: nil
  },
  'nad5' => {
    size: nil
  },
  'nad6' => {
    size: nil
  },
  'cob' => {
    size: nil
  }
}

outgroup = 'NC_025320_1'
bootstraps = 1000
rows = {
  'KC990830_1' => [], # agilis
  'NC_011607_1' => [], # podarcis
  'NC_018777_1' => [], # takydromus
  'NC_025320_1' => [], # eremias
  'NC_026867_1' => [], # zootoca
  'boehme_lacerta_viridis' => [],
  'rohit_adriatic_lineage' => [],
  'rohit_lacerta_bilineata' => [],
  'rohit_lacerta_viridis' => [],
  'S1_S1' => [],
  'S2_S22' => [],
  'S3_S16' => [],
  'S4_S23' => [],
  'S5_S5' => [],
  'S6_S6' => [],
  'S7_S13' => [],
  'S8_S14' => [],
  'S9_S9' => [],
  'S10_S10' => [],
  'S11_S11' => [],
  'S12_S21' => [],
  'S17_S17' => [],
  'S18_S18' => [],
  'S19_S19' => [],
  'S24_S24' => []
}

TAXA_NAMES = {
'S1_S1' => 'G841_384_Lb_FR',
'S2_S22' => 'G664_157_Lv_CZ',
'S3_S16' => 'G2610__Lv_BG',
'S4_S23' => 'G593_72_Lb_ES',
'S5_S5' => 'G2029_89_Lv_BG',
'S6_S6' => 'G2167_252_Lv_BG',
'S7_S13' => '9003__Lv_TR',
'S8_S14' => '9004__Lv_TR',
'S9_S9' => 'G380_204_AL_ME',
'S10_S10' => 'G707_216_Lv_RS',
'S11_S11' => 'G899_449_Lv_MD',
'S12_S21' => 'G2697__Lb_IT',
'S17_S17' => 'G2706__Lb_FR',
'S18_S18' => 'G1993_53_Lv_UA',
'S19_S19' => 'G2705__Lb_FR',
'S24_S24' => 'G376_397_AL_GR',
'boehme_lacerta_viridis' => 'AM176577_1_Lv_AT',
'rohit_lacerta_viridis' => 'Lv_HU',
'rohit_adriatic_lineage' => 'AL_SI',
'rohit_lacerta_bilineata' => 'KT722705__Lb_FR',
'KC990830_1' => 'KC990830__La',
'NC_010972_2' => 'NC010972_2_Ac',
'NC_011607_1' => 'NC_011607_1_Pm',
'NC_018777_1' => 'NC_018777_1_Tw',
'NC_025320_1' => 'NC_025320_1_Ev',
'NC_026867_1' => 'NC_026867_1_Zv'
}

partition = []
partition_ng = []
partition_nexus = ['begin mrbayes;', '  set autoclose=yes;', "  outgroup #{outgroup};"]
start = 1
gen_hash.keys.each do |gen|
  al = Bio::Alignment::MultiFastaFormat.new(IO.read("./#{base_dir}_genes/#{gen}_maffted.fas"))
  size = nil
  al.entries.each do |fasta|
    next unless rows.keys.include?(fasta.definition)
    gen_hash[gen][:size] = fasta.seq.size
    rows[fasta.definition] << fasta.seq
  end
  partition << "DNA, #{gen} = #{start}-#{start + gen_hash[gen][:size] -1}"
  partition_ng << "HKY+G, #{gen} = #{start}-#{start + gen_hash[gen][:size] -1}"
  partition_nexus << "  charset #{gen} = #{start}-#{start + gen_hash[gen][:size] -1};"
  start += gen_hash[gen][:size]
end
# write raxml patition
IO.write("./#{base_dir}_genes/raxml_partition.txt", partition.join("\n"))
# write default raxml_ng partition (all with HKY+G, must be adjusted)
IO.write("./#{base_dir}_genes/raxml_ng_partition.txt", partition_ng.join("\n"))
# write nexus partition and block
partition_nexus << "  partition genes = #{gen_hash.keys.size}: #{gen_hash.keys.join(', ')};"
partition_nexus << '  set partition = genes;'
partition_nexus << '  lset applyto=(all) rates=gamma;'
partition_nexus << '  prset applyto=(all) ratepr=variable;'
partition_nexus << '  mcmcp samplefreq=100 printfreq=500;'
partition_nexus << '  mcmcp nruns=2 stoprule=YES burninfrac=.25;'
partition_nexus << '  mcmcp stopval=0.01 minpartfreq=0.05;'
partition_nexus << '  mcmcp mcmcdiagn=YES diagnfreq=5000;'
partition_nexus << '  mcmc;'
partition_nexus << '  sump;'
partition_nexus << '  sumt Conformat=Simple;'
# partition_nexus << '  mcmc ngen=1000000 samplefreq=500 printfreq=500 diagnfreq=5000;'
partition_nexus << 'end;'
IO.write("./#{base_dir}_genes/mrbayes.partition.txt", partition_nexus.join("\n"))

complete = ''
rows.each_pair do |k,v|
  complete += ">#{k}\n"
  complete += v.join
  complete += "\n"
end
# write fasta file
IO.write("./#{base_dir}_genes/raxml_complete.fas", complete)
# write nexus file
Dir.chdir("./#{base_dir}_genes") do
  `../fasta2nex.r raxml_complete.fas mrbayes.alignment.nexus`
  `cat mrbayes.alignment.nexus mrbayes.partition.txt > reconstruct/mrbayes_complete_#{TAXA_NAMES[outgroup]}.nexus`
end

def convert_tree_names(tree, stringify)
  org = IO.read(tree)
  TAXA_NAMES.each_pair do |k, v|
    if org.include?(k)
      if stringify
        org.gsub!(k, "'#{v}'")
      else
        org.gsub!(k, v)
      end
    else
      puts "missing #{v}"
    end
  end
  IO.write(tree, org)
end

def convert_nexus_to_nwk(nexus_tree)
  nexus = Bio::Nexus.new( IO.read( nexus_tree ) )
  trees_block = nexus.get_trees_blocks[ 0 ]
  tree_first = trees_block.get_tree( 0 )

  taxa = {}
  if nexus.get_taxa_blocks[ 0 ]
    taxa_block = nexus.get_taxa_blocks[ 0 ]
    i = -1
    puts "TAXA Block"
    taxa = Hash[[nil, taxa_block.get_taxa].flatten.collect{|e| [(i += 1).to_s, e]}]
  else
    puts "TRANSLATE Block"
    taxa = Hash[content.match(/Translate\s*([^;]*);/i)[1].split(',').map(&:strip).collect{|e| e.split(' ')}]
  end

  tree_first.each_node{ |n|
    next if n.name.to_i == 0
    n.name = taxa[n.name]
  }
  tree = tree_first.output_newick
  if tree.match(/^\(\s*\[\&U\],\s*/i)
    tree = tree.gsub(/^\(\s*\[\&U\],\s*/i, '').gsub(/\s*\);/,';')
  end
  nwk_tree = nexus_tree.gsub('mrbayes_complete_', 'final_tree_mrbayes_')
  nwk_tree.gsub!('.nexus.con.tre', '.nwk')
  IO.write(nwk_tree, tree)
end

FileUtils.mkdir_p "./#{base_dir}_genes/reconstruct"
Dir.chdir("./#{base_dir}_genes/reconstruct") do
  # --- RAxML (standard) ---
  # run binary
  system "/bio-bin/standard-RAxML/raxmlHPC-SSE3 -f a -s ../raxml_complete.fas -n lacerta_13_gene_#{bootstraps} -m GTRGAMMA -x 1234 -##{bootstraps} -p 1234 -q ../raxml_partition.txt -o #{outgroup}"
  system "cp RAxML_bipartitions.lacerta_13_gene_#{bootstraps} final_tree_raxml_#{bootstraps}_gtrgamma_#{TAXA_NAMES[outgroup]}.nwk"
  # rename tree
  convert_tree_names("final_tree_raxml_#{bootstraps}_gtrgamma_#{TAXA_NAMES[outgroup]}.nwk", true)
  # --- MrBayes ---
  # run binary
  system "/bio-bin/MrBayes/mb mrbayes_complete_#{TAXA_NAMES[outgroup]}.nexus"
  # convert nexus tree to wnk
  convert_nexus_to_nwk("mrbayes_complete_#{TAXA_NAMES[outgroup]}.nexus.con.tre")
  # rename tree
  convert_tree_names("final_tree_mrbayes_#{TAXA_NAMES[outgroup]}.nwk", false)
end
