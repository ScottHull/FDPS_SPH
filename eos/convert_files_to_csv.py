import os
import csv

fname1 = 'granite.rho_u.txt'
fname2 = 'granite.table.txt'
fnames = [fname1, fname2]

def returnReader(fname, delimiter):
    with open(fname, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        return list(reader)

def removeOutput(oname):
    if oname in os.listdir(os.getcwd()):
        os.remove(oname)


for fname in fnames:
    oname = fname.replace('.txt', '.csv')
    removeOutput(oname=oname)
    with open(oname, 'w') as outfile:
        reader = returnReader(fname=fname, delimiter=' ')
        for row in reader:
            if len(row) == 0:
                pass
            else:
                rowstr = ",".join(row).replace(",,", ',')
                try:
                    if rowstr[0] == ',':
                        rowstr = rowstr[1:]
                        outfile.write(rowstr + '\n')
                    else:
                        outfile.write(rowstr + '\n')
                except:
                    outfile.write(rowstr + '\n')
    outfile.close()