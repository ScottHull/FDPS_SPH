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


oname = fname1.replace('.txt', '.csv')
removeOutput(oname=oname)
with open(oname, 'w') as outfile:
    reader = returnReader(fname=fname1, delimiter='\t')
    for row in reader:
        rowstr = ",".join(row)
        outfile.write(rowstr + '\n')
outfile.close()

oname = fname2.replace('.txt', '.csv')
removeOutput(oname=oname)
with open(oname, 'w') as outfile:
    reader = returnReader(fname=fname2, delimiter=' ')
    for row in reader:
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
