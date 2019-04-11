#word to phones
import sys

try:
        filename=sys.argv[1]
except IndexError:
        print("usage: python generate_lxicon.py <path_to_vocabulary_file>")
	print("e.g python generate_lxicon.py words.txt")
	print("Generate the lexicon.txt in the present directory")
        sys.exit(1)

from g2p_en import g2p
import os

filepath=os.path.dirname("filename")
print(filename)
myfile=open("lexicon.txt", "w") 
with open(filename) as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content] 
for x in content:
	myfile.write(x)
	for i in g2p(x):
		myfile.write(" ")
		myfile.write(str(i))
	myfile.write("\n")
myfile.close()
