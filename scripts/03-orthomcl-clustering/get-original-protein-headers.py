import sys

proteome_dict = {}
with open(sys.argv[2], 'r') as f:
    for line in f:
        line = line.strip('\n')
        line = line.split('\t')
        proteome_dict[line[0]] = line[1]

with open(sys.argv[1], 'r') as f:
    og_counter = len(f.readlines())

output = open(sys.argv[3], 'w')

with open(sys.argv[1], 'r') as f:
    for i, line in enumerate(f):
        line = line.strip('\n')
        line = line.split(' ')
        n=len(list(str(og_counter))) + 2
        new_line = ['OG' + str(i).zfill(n) + ':']
        for id in line[1:]:
            new_line.append(proteome_dict[id])
        new_line = ' '.join(new_line) + '\n'
        output.write(new_line)

