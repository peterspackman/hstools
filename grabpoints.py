import sys

vertices = 0
linesout = []

with open(sys.argv[1]) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if vertices > 0:
            linesout.append(line)
            vertices = vertices - 1
        if(line.strip().startswith('begin vertices_recon')):
            vertices = int(line.split()[2])
            print 'Found vertices: {0}'.format(vertices)
        if vertices < 0:
            break


with open(sys.argv[2], 'w') as f:
    for line in linesout:
        f.write(line)
