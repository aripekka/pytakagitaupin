#Extracts the susceptibility data from xfh.inf

xfh_file=open('xfh.inf')

energy=[]
chi=[]

for line in xfh_file:
    lines = line.split()

    if len(lines) > 1:

        #get energy in keV
        if lines[0] == 'at' and lines[1] == 'energy':
            energy.append(float(lines[3])/1000)

        #get chi_0
        elif lines[0] == 'Psi_0':
            aux = lines[2]
            aux = aux[1:-1].split(',')

            chi.append([])
            chi[-1].append(-float(aux[0]))
            chi[-1].append(-float(aux[1]))

        #get chi_h
        elif lines[0] == 'Psi_H':
            aux = lines[2]
            aux = aux[1:-1].split(',')

            chi[-1].append(-float(aux[0]))
            chi[-1].append(-float(aux[1]))

        #get chi_hbar
        elif lines[0] == 'Psi_HBar':
            aux = lines[2]
            aux = aux[1:-1].split(',')

            chi[-1].append(-float(aux[0]))
            chi[-1].append(-float(aux[1]))

xfh_file.close()

#write the chitable file

chi_file = open('chitable.dat','w')

for i in range(len(energy)):
    line = str(energy[i])
    for chival in chi[i]:
        line = line + ' ' + str(chival)
    chi_file.write(line + '\n')

chi_file.close()

