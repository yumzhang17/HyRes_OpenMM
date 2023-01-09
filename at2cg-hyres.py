#!/usr/bin/python

import sys

if len(sys.argv) != 3:
   print "Usage: at2cg.py inp.pdb out.pdb"
   print "       in the inp.pdb, it's atomistic structure"
   print "                       all HSD, HSE, HSP has to be renamed to HIS"
   print "                       atom name and resname follows GBSW model"
   print "       in the out.pdb, it's IDPCG structure"
   print "       in the mapping process, geometric COM is computed"
   quit()

inp=sys.argv[1]
out=sys.argv[2]

f1=open(inp,'r')
f2=open(out,'w')

# read input pdb file
data={} 
tmp={}
natom=0
nres=1
pre=1
iatom=0
for line in f1:
    if line.startswith("REMARK"):
       pass
    elif line.startswith("ATOM"):
         natom=natom+1
         if line.split()[4] == str(pre):
            iatom=iatom+1
            tmp[iatom]=line.split()
         elif line.split()[4] != str(pre):
            data[nres]=tmp
            nres=nres+1
            pre=line.split()[4]
            iatom=1
            tmp={}
            tmp[iatom]=line.split()
    else:
       pass
    data[nres]=tmp
print "There are",natom,"atoms /",nres,"residues"

# mapping rules
reslist=['gly', 'ala', 'val', 'leu', 'ile', 'met', 'asn', 'asp', 'gln', 'glu', 'cys', 'ser', 'thr', 'pro', 'lys', 'arg', 'his', 'phe', 'tyr', 'trp']
single=['ala', 'val', 'leu', 'ile', 'met', 'asn', 'asp', 'gln', 'glu', 'cys', 'ser', 'thr', 'pro']
def maprule(resname):
    nsc=0
    sc1=[]
    sc2=[]
    sc3=[]
    sc4=[]
    sc5=[]
    bb=[]
    nter=[]
    cter=[]
    if resname in reslist:
       nter=['CAY', 'HY1', 'HY2', 'HY3', 'CY', 'OY']
       cter=['NT', 'HNT', 'CAT', 'HT1', 'HT2', 'HT3']
    if resname in reslist and resname != 'pro' and resname != 'gly':
       bb=['CA', 'HA', 'C', 'O', 'N', 'HN']
    elif resname == 'gly':
       bb=['CA', 'HA1', 'HA2', 'C', 'O', 'N', 'HN']
    elif resname == 'pro':
       bb=['CA', 'HA', 'C', 'O', 'N']
    if resname == 'lys':
       nsc=2
       sc1=['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2']
       sc2=['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2', 'HZ3']
    elif resname == 'arg':
       nsc=2
       sc1=['CB', 'HB1', 'HB2', 'CG', 'HG1', 'HG2', 'CD', 'HD1', 'HD2']
       sc2=['NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22']
    elif resname == 'his':
       nsc=3
       sc1=['CB', 'HB1', 'HB2', 'CG']
       sc2=['CD2', 'HD2', 'NE2']
       sc3=['ND1', 'HD1', 'CE1', 'HE1']
    elif resname == 'phe':
       nsc=3
       sc1=['CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1']
       sc2=['CD2', 'HD2', 'CE2', 'HE2']
       sc3=['CE1', 'HE1', 'CZ', 'HZ']
    elif resname == 'tyr':
       nsc=3
       sc1=['CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1']
       sc2=['CD2', 'HD2', 'CE2', 'HE2']
       sc3=['CE1', 'HE1', 'CZ', 'OH', 'HH']
    elif resname == 'trp':
       nsc=5
       sc1=['CB', 'HB1', 'HB2', 'CG']
       sc2=['CD1', 'HD1', 'NE1', 'HE1']
       sc3=['CD2', 'CE2']
       sc4=['CZ2', 'HZ2', 'CH2', 'HH2']
       sc5=['CE3', 'HE3', 'CZ3', 'HZ3']
    elif resname in single:
       nsc=1
    elif resname == 'gly':
       nsc=0
    return nsc,bb,nter,cter,sc1,sc2,sc3,sc4,sc5

# output in pdb-format
def printcg(atom):
    f2.write("%4s   %4d %4s %4s  %3d     %7.3f %7.3f %7.3f  %4.2f %5.2f      %4s\n" % (atom[0],int(atom[1]),atom[2],atom[3],int(atom[4]),float(atom[5]),float(atom[6]),float(atom[7]),float(atom[8]),float(atom[9]),atom[10]))

# convert at to cg
bbcg1=['CA', 'N', 'HN']
bbcg2=['C', 'O']
ntercg=['CAY', 'CY', 'OY']
ctercg=['NT', 'HNT', 'CAT']
inx=0
for ires in range(1,nres+1,1):   # loop over all residues
    iresname=data[ires][1][3].lower()  # residue number
    iresatnum=len(data[ires].keys())   # number of atoms in the residue
    if iresname in reslist:
       # mapping rules
       (nsc,bb,nter,cter,sc1,sc2,sc3,sc4,sc5)=maprule(str(iresname))
       if iresname in single:
          sc1=[item[2] for item in data[ires].values() if item[2] not in bb if item[2] not in nter if item[2] not in cter]
       # initializing
       coor={}
       num={}
       for i in range(1,6,1):
           coor[i]=[0,0,0]
           num[i]=0
       # compute com for each sc bead
       for j in range(1,iresatnum+1,1): # loop over all atoms in the residue
           if data[ires][j][2] in ntercg:
              data[ires][j][3]=data[ires][j][3]+'_'
           elif data[ires][j][2] in ctercg:
              data[ires][j][3]=data[ires][j][3]+'_'
           elif data[ires][j][2] in bbcg1:
              inx=inx+1
              if data[ires][j][2] == 'HN':
                 data[ires][j][2]='H'
              data[ires][j][1]=inx
              data[ires][j][3]=data[ires][j][3]+'_'
              printcg(data[ires][j])
           elif data[ires][j][2] in sc1:
              num[1]=num[1]+1
              coor[1][0]=coor[1][0]+float(data[ires][j][5]) # x
              coor[1][1]=coor[1][1]+float(data[ires][j][6]) # y
              coor[1][2]=coor[1][2]+float(data[ires][j][7]) # z
           elif data[ires][j][2] in sc2:
              num[2]=num[2]+1
              coor[2][0]=coor[2][0]+float(data[ires][j][5]) # x
              coor[2][1]=coor[2][1]+float(data[ires][j][6]) # y
              coor[2][2]=coor[2][2]+float(data[ires][j][7]) # z
           elif data[ires][j][2] in sc3:
              num[3]=num[3]+1
              coor[3][0]=coor[3][0]+float(data[ires][j][5]) # x
              coor[3][1]=coor[3][1]+float(data[ires][j][6]) # y
              coor[3][2]=coor[3][2]+float(data[ires][j][7]) # z
           elif data[ires][j][2] in sc4:
              num[4]=num[4]+1
              coor[4][0]=coor[4][0]+float(data[ires][j][5]) # x
              coor[4][1]=coor[4][1]+float(data[ires][j][6]) # y
              coor[4][2]=coor[4][2]+float(data[ires][j][7]) # z
           elif data[ires][j][2] in sc5:
              num[5]=num[5]+1
              coor[5][0]=coor[5][0]+float(data[ires][j][5]) # x
              coor[5][1]=coor[5][1]+float(data[ires][j][6]) # y
              coor[5][2]=coor[5][2]+float(data[ires][j][7]) # z
       for i in range(1,nsc+1,1):
           inx=inx+1
           if i == 1:
              name='CB'
           elif i == 2:
              name='CC'
           elif i == 3:
              name='CD'
           elif i == 4:
              name='CE'
           elif i == 5:
              name='CF'
           coor[i][0]=coor[i][0]/num[i]
           coor[i][1]=coor[i][1]/num[i]
           coor[i][2]=coor[i][2]/num[i]
           tmp=[data[ires][1][0],inx,name,data[ires][1][3],data[ires][1][4],coor[i][0],coor[i][1],coor[i][2],data[ires][1][8],data[ires][1][9],data[ires][1][10]]
           printcg(tmp)
       # this is to make sure cg atoms are in correct order
       for j in range(1,iresatnum+1,1): # loop over all atoms in the residue
           if data[ires][j][2] in bbcg2:
              inx=inx+1
              data[ires][j][1]=inx
              data[ires][j][3]=data[ires][j][3]+'_'
              printcg(data[ires][j])
    else:
       print iresname,"is not recognized" 
       quit()
          
f2.write("%3s\n" % ("END"))
quit()
