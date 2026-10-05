import sys,os


bed=sys.argv[1]
lib=sys.argv[2]

library={}
new_deletion={}
for i in open(lib):
    i=i.rstrip()
    if i not in library:
        library[i]="Y"

for j in open(bed):
    j=j.rstrip()
    if j not in library:
        library[j]="Y"
        new_deletion[j]="Y"

#update library and bed file
os.system('mv '+lib+" "+lib+".prev")
os.system('mv '+bed+" "+bed+".prev")

out_lib=open(lib,'w')
if len(library) != 0 :
    for lib in library:
            out_lib.write(lib+'\n')
out_lib.close()

out_bed=open(bed,'w')
if len(new_deletion) != 0 :
    for de in new_deletion:
            out_bed.write(de+'\n')
out_bed.close()
