import regina
import sys

t = regina.Triangulation3(sys.argv[1])
t.orient()

h = eval(sys.argv[2])

expanded = (4 * t.countTetrahedra()) * [ 0 ]

for f, c in zip(t.faces(2), h):
    for i, e in enumerate(f.embeddings()):
        expanded[4 * e.tetrahedron().index() + e.vertices()[3]] += (-1) ** i * c

print(expanded)
