import sys

with open(sys.argv[1], 'r') as in_f, open(sys.argv[2], 'w') as out_f:
  line_counter = 0

  for line in in_f:
    if line_counter == 0:
      toks = line.split()

      idx1 = toks.index('HGDP01006')
      del toks[idx1]
      del toks[idx1]
      idx2 = toks.index('HGDP01051')
      del toks[idx2]
      del toks[idx2]

      inds = toks[2:]
      c_inds = ','.join(inds)
      out_line = f'#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,{c_inds}\n'
      out_f.write(out_line)
    else:
      toks = line.split()

      del toks[idx1]
      del toks[idx1]
      del toks[idx2]
      del toks[idx2]

      site = toks[0]
      GTs  = toks[2:]
      c_GTs = ','.join(GTs)
      out_line = f'6,{site},.,{GTs[0]},-,.,PASS,TARGT,GT,{c_GTs}\n'
      out_f.write(out_line)
    line_counter += 1
