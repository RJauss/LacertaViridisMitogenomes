#!/usr/bin/env ruby
# frozen_string_literal: true

require 'bio'

genes = /(atp6|atp8|cob|cox1|cox2|cox3|nad1|nad2|nad3|nad4l|nad4|nad5_0|nad5_1|nad5|nad6)/

base_dir = 'MITOS_2019_08_09'

system "rm -rf ./#{base_dir}_genes"

Dir["#{base_dir}/*"].each do |dir|
  file = Dir["#{dir}/*.fas"].first
  next if file.nil?

  al = Bio::Alignment::MultiFastaFormat.new(IO.read(file))
  al.entries.each do |fasta|
    header = fasta.definition.split('; ')
    next unless header.last.match genes

    gen = header.last.match(genes)[1]
    fasta.definition = File.basename(dir)
    system 'mkdir', '-p', "./#{base_dir}_genes/#{gen}"
    IO.write("./#{base_dir}_genes/#{gen}/#{fasta.definition}.fas", fasta.to_s)
    IO.write("./#{base_dir}_genes/nad5/#{fasta.definition}.fas", fasta.to_s) if gen == 'nad5_0'
  end
end

Dir["#{base_dir}_genes/*"].each do |dir|
  gen = File.basename(dir)
  system "cat ./#{base_dir}_genes/#{gen}/*.fas >./#{base_dir}_genes/#{gen}_unaligned.fas"
  system "mafft --localpair ./#{base_dir}_genes/#{gen}_unaligned.fas >./#{base_dir}_genes/#{gen}_maffted.fas"
end
