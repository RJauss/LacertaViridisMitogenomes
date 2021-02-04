#!/usr/bin/env ruby
# frozen_string_literal: true

require 'bio'
require 'pp'
require 'descriptive_statistics'

base_dir = 'MITOS_2019_08_09'

CURRENT_NUMBER = 14

gen_hash = {
  'nad1' => {},
  'nad2' => {},
  'cox1' => {},
  'cox2' => {},
  'atp8' => {},
  'atp6' => {},
  'cox3' => {},
  'nad3' => {},
  'nad4l' => {},
  'nad4' => {},
  'nad5' => {},
  'nad6' => {},
  'cob' => {}
}

TAXA_NAMES = {
  'boehme_lacerta_viridis' => 'AM176577_1_Lv_AT',
  'rohit_lacerta_viridis' => 'Lv_HU',
  'S2_S22' => 'G664_157_Lv_CZ',
  'S3_S16' => 'G2610__Lv_BG',
  'S5_S5' => 'G2029_89_Lv_BG',
  'S6_S6' => 'G2167_252_Lv_BG',
  'S10_S10' => 'G707_216_Lv_RS',
  'S11_S11' => 'G899_449_Lv_MD',
  'S18_S18' => 'G1993_53_Lv_UA',
  'S7_S13' => '9003__Lv_TR',
  'S8_S14' => '9004__Lv_TR',
  'S9_S9' => 'G380_204_AL_ME',
  'S24_S24' => 'G376_397_AL_GR',
  'S12_S21' => 'G2697__Lb_IT',
  'rohit_adriatic_lineage' => 'AL_SI',
  'rohit_lacerta_bilineata' => 'KT722705__Lb_FR',
  'S1_S1' => 'G841_384_Lb_FR',
  'S4_S23' => 'G593_72_Lb_ES',
  'S17_S17' => 'G2706__Lb_FR',
  'S19_S19' => 'G2705__Lb_FR',
}

MOTIFS = {
  'CSB-1' => 'ctatatggtattattgtcttaatgcttggtagacatat',
  'CSB-2' => 'caaacccccctacccccc',
  'CSB-3' => 'tcgccaaacccctaaaacga',
  'R1.*' => 'caccactttctcactttttccaaggcctctggtt',
  'R1.a' => 'catcactttctcactttttccaaggcctctggtt',
  'R1.b' => 'caccactttctcaccttttccaaggcctctggtt',
  'R1.c' => 'caccactttatcaccttttccaaggcctctggtt',
  'R2' => 'atttaccccatgaa',
  'TAS1' => 'actattatgtatatagtgcattaa',
  'TAS2.*' => 'catacattaa',
  'TAS2.a' => 'tatgcattaa',
  'TAS2.b' => 'tatacattaa',
  # 'OH' => 'tatatactattatgtatatagtgcattaactgatttaccccatgaaaaataatatgtacattatcctattaataagacataatacatatatgtataatcatacattaatatatttaccccatgaatatccttgtatatactagtatctatatattttacataatacataactgtttaaagtacatagcattataaaccatacgactatt',
  'Lv_Insertion' => 'aacaaactttaa',
  'rep35bp.*' => 'cacctgccgcttaaaagcggcttttttgcctccta',
  'rep35bp.1' => 'cacccgccgcttaaaaagcggcttttttgcctccta', # created at pos 0
  'rep35bp.2' => 'cacccgccgctaactagcggcttttttgcctccta', # created at pos 0
  'rep35bp.3' => 'cacccgccgcttaaaagcggcttttttgcctccta', # created at pos 0
  'rep35bp.a' => 'cacctgccgcttaaaagcggctttttgcctccta', # created at pos 0
  'rep35bp.b' => 'cacccgccgctttaaagcggcttttttgcctccta', # created at pos 0
  'rep35bp.c' => 'cacccgccgctaacgcggcttttttgcctccta', # created at pos 0
  'rep35bp.d' => 'cacccgccgctaactagcggcttttttgcctcctg', # created at pos 0
  'rep35bp.e' => 'cacccgccgcttagaaagcggcttttttgcctccta', # created at pos 0
  'rep35bp.f' => 'cacccgccgcttaaaaagcggcttttttgcctcctg', # created at pos 0
  'rep35bp.g' => 'cacccgccgcttagtcttttttgcaggtgc', # 9003_Lv_TR
  'rep35bp.h' => 'cacctgccgcttagaagcggcttttttgcctccta', # created for pos 4
  'rep35bp.i' => 'cacccgccgcttagaagcggcttttttgcctccta', # created for pos 4
  'rep35bp.j' => 'cacccgccgcttaaaaagcggcttttttgcctcctccta', # created for pos 4 (G380_204_AL_ME)
  'rep35bp.k' => 'cacctgccgcttagaagcggcttttctgcctccta', # created for pos 5 (G2610__Lv_BG)
  'rep35bp.l' => 'ctcccgccgctttaaagcggcttttttgcctccta', # created for pos 5 (9003__Lv_TR)
  'rep35bp.m' => 'cgcccgccgcttaaaaagcggcttttttgcctccta', # created for pos 6 (G380_204_AL_ME)
  'rep35bp.n' => 'tacccgccgcttaaaaagcggcttttttgcctccta', # created for pos 5 (G380_204_AL_GR) ATTENTION: leading T is second last from prev 35er motif (2 gap in a row)
  'rep35bp.o' => 'cacccgccgcttagtactttttaaagcggcttttttgcctccta', # created for pos 5 (AL_SI)
  'rep35bp.p' => 'caccccccggtaactagcggtttttttgccctccta', # created for pos 5 (G841_384_Lb_FR)
  'rep35bp.q' => 'cacccgccgctaactagcggcttttttgccccctc', # created for pos 5 (G2706__Lb_FR)
  'rep35bp.r' => 'ccccccccgctaactagcggcttttttgcctccta', # created for pos 6 (G2706__Lb_FR)
  'rep35bp.s' => 'cacccgccgctatcgcggcttttttgcctcct', # created for pos 6 (G593_72_Lb_ES)
  'rep35bp.t' => 'cacctgccgcttaaaagcggcttttttgccttctg', # created for pos 6 (G2167_252_Lv_BG)
  'rep35bp.u' => 'cccccgccgctaactagcggcttttttgcctccta', # created for pos 6 (G841_384_Lb_FR)
  'rep35bp.v' => 'cacccgccccaaagttgccgtttttttcccccccac', # created for pos 5 (G2705__Lb_FR)
  'rep35bp.w' => 'ccccccccgttagttggcgttttttttccctcctc', # created for pos 6 (G2705__Lb_FR)
  'rep35bp.x' => 'cccccgccgctaacgagcgtcttttttgcctcctt', # created for pos 7 (G2705__Lb_FR)
  'rep35bp.y' => 'cacccgccgcttaaaaagcggcttttttgcctccg', # created for pos 8 (9003__Lv_TR)
  'rep35bp.z' => 'caccnnnnnnnnnnnnnnnnnnnnnnngcctccta' # created at pos 0
}
rep35bp =  MOTIFS.collect do |k,v|
  [">#{k}", v] if k.include?('rep35bp')
end
IO.write("dloop.alignment.#{CURRENT_NUMBER}.fas", rep35bp.compact.join("\n"))
`mafft --localpair dloop.alignment.#{CURRENT_NUMBER}.fas > dloop.alignment.#{CURRENT_NUMBER}.maffted.fas`

alignment = Bio::Alignment::MultiFastaFormat.new(IO.read("dloop.alignment.#{CURRENT_NUMBER}.maffted.fas"))
rep35bp_info = {}
ref = alignment.entries.first.seq
alignment.entries.each do |am|
  rep35bp_info[am.definition.gsub('rep35bp.', '')] = {
    'motif' => am.seq,
    'diff' => am.seq.chars.zip(ref.chars).select { |a,b| a != b }.count,
    'perc' => (am.seq.chars.zip(ref.chars).select { |a,b| a != b }.count * 100.0 / 35.0).round(2),
    'found' => 0
  }
end

animals = {}
fasta = []
TAXA_NAMES.keys.each do |k|
  arr = IO.readlines(Dir["#{base_dir}/#{k}/*.bed"].first, chomp: true).collect{|l| l.split("\t")[1...-2]}
  arr = arr.reject{|e| !e[2].include?('trnP') && !e[2].include?('trnF') }.collect{|e| [e[0].to_i, e[1].to_i, e[2]]}
  arr = arr.sort{|r,l| l.last <=> r.last}
  start, stop = [arr.first[1], arr.last[1]]
  seq = Bio::Alignment::MultiFastaFormat.new(IO.read(Dir["#{base_dir}/#{k}/*.fasta"].first)).alignment.entries.first.seq
  dloop = ''
  if (start < stop)
    dloop = seq[start...stop]
  else
    dloop = seq[start...-1]
    dloop+= seq[0...stop]
  end
  animals[TAXA_NAMES[k]] = {
    'contig' => seq.size,
    'start' => start,
    'stop' => stop,
    'dloop' => dloop,
    'motifs' => {},
    'rep35bp' => {}
  }
  MOTIFS.each_pair do |name, motif|
    animals[TAXA_NAMES[k]]['motifs'][name] = []
    ofs = 0
    while(index = animals[TAXA_NAMES[k]]['dloop'].index(motif, ofs)) do
      animals[TAXA_NAMES[k]]['motifs'][name] << index
      rep35bp_info[name.gsub('rep35bp.', '')]['found'] += 1 if name.include?('rep35bp.')
      animals[TAXA_NAMES[k]]['rep35bp'][index] = name.gsub('rep35bp.', '') if name.include?('rep35bp.')
      ofs = index + 1
    end
  end
  # puts "#{TAXA_NAMES[k]}\tdloop: #{dloop.size} pb"
  fasta << "> #{TAXA_NAMES[k]}"
  fasta << dloop
end

File.open("dloop.35er-motif.#{CURRENT_NUMBER}.tsv", 'w') do |f|
  f.puts "Taxa\tstart\tstop\tA\tT\tG\tC\tA+T\tG+C"
  nucs_total = [
    [],
    [],
    [],
    [],
    [],
    []
  ]
  animals.each_pair do |taxa, values|
    s = values['rep35bp'].keys.sort.first
    e = values['rep35bp'].keys.sort.last
    m = values['rep35bp'][e]
    e += MOTIFS["rep35bp.#{m}"].size
    repl = values['dloop'][s...e]
    nucs = [
      repl.count('aA'),
      repl.count('tT'),
      repl.count('gG'),
      repl.count('cC'),
    ]
    total = nucs.sum
    nucs.map! { |e| (e * 100.0 / total).round(3) }
    nucs << (nucs[0] + nucs[1]).round(3)
    nucs << (nucs[2] + nucs[3]).round(3)
    nucs.each_with_index do |p, i|
      nucs_total[i] << p unless taxa.include?('G2029_89_Lv_BG') # excluded because of N
    end
    f.print taxa
    f.print "\t"
    f.print s
    f.print "\t"
    f.print e
    f.print "\t"
    f.print nucs.join("\t")
    f.puts
  end

  tl = ['Median', '-', '-']
  nucs_total.each do |perc|
    stats = [
      perc.min,
      perc.max,
      perc.percentile(50).round(3),
      perc.percentile(25).round(3),
      perc.percentile(75).round(3),
      (perc.percentile(75) - perc.percentile(25)).round(3),
      (perc.percentile(25) - 1.5 * (perc.percentile(75) - perc.percentile(25))).round(3),
      (perc.percentile(75) + 1.5 * (perc.percentile(75) - perc.percentile(25))).round(3),
      (perc.percentile(25) - 3.0 * (perc.percentile(75) - perc.percentile(25))).round(3),
      (perc.percentile(75) + 3.0 * (perc.percentile(75) - perc.percentile(25))).round(3)
    ]
    tl << stats[2]
  end
  f.puts tl.join("\t")
end


content = ''
content = "Taxa\tcontig\td-loop\t" + MOTIFS.keys.join("\t") + "\trep.count\trep.pos\trep.types\n"
animals.each_pair do |taxa, values|
  content += "#{taxa}\t#{values['contig']}\t#{values['dloop'].size}\t"
  totals = []
  values['motifs'].each_pair do |motif, positions|
    content += positions.empty? ? '-' : positions.join('|')
    content += "\t"
    totals += positions if motif.include?('rep35bp')
  end
  content += "#{totals.uniq.size}"
  content += "\t"
  content += totals.uniq.empty? ? '-' : totals.uniq.sort.join('|')
  content += "\t"
  content += values['rep35bp'].sort.collect{|e| e.last}.join
  content += "\n"
end
IO.write("dloop.table.#{CURRENT_NUMBER}.tsv", content)

File.open("dloop.info.#{CURRENT_NUMBER}.tsv", 'w') do |f|
  f.puts "type\tmotif\tdiff\tpercent\tfound"
  rep35bp_info.each_pair do |k, v|
    f.print k
    f.print "\t"
    f.print v['motif']
    f.print "\t"
    f.print v['diff']
    f.print "\t"
    f.print v['perc']
    f.print "\t"
    f.print v['found']
    f.puts
  end
end
