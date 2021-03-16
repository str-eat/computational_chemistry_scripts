lens = ""
angs = ""
dihs = ""

outfile = input("Enter output filename: ")
infile = input("Enter input filename: ")

with open(outfile, "w+") as w:
	with open(infile, "r+") as r:
		for line in r:
			try:
				if not line.split()[-1][0].isnumeric():
					raise Exception
				elif line.split()[0][0:3] == 'dih':
					i = 4
				elif line.split()[0][2].isnumeric():
					i = 2
				else:
					i = 3

				if i == 2:
					print(i)
					lens = lens + "{}    {:0.1f}\n".format(line.split()[0], float(line.split()[-1]))
				elif i == 3:
					print(i)
					angs = angs + "{}    {:0.0f}.0\n".format(line.split()[0], float(line.split()[-1]))
				elif i == 4:
					print(i)
					dihs = dihs + "{}    {:0.0f}.0\n".format(line.split()[0], float(line.split()[-1]))
			except:
				try:
					if line.split()[0][0:3] == 'dih':
						dihs = dihs + "{}    {:0.0f}.0\n".format(line.split()[0], float(line.split()[-1]))
					else:
						print(line)
						w.write(line)
				except:
					print(line)
					w.write(line)
	w.write(lens + angs + dihs)
