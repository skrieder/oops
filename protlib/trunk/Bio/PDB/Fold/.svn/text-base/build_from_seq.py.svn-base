#!/usr/bin/python
"""This is a module to take a protein sequence as stored in a suitable file type (suggested FASTA) and generate a three dimensional conformation"""
from Bio import SeqIO,PDB
import StringIO
import sys,os,math

NHLength=1 #angstrom

try:
    import numpy
except ImportError:
    def CUMSUM(a):
        l = [a[0]]
        for i in xrange(1, len(a)):
            l.append(a[i] + l[-1])
        return l
    def ZEROS(n, dummy):
        return [0.0]*n
    def ARRAY(v, dummy):
        return list(v)
else:
    CUMSUM = numpy.add.accumulate
    ZEROS = numpy.zeros
    ARRAY = numpy.array

class Atom:
    def __init__(self, name, p1, p2, p3, d1, d2, d3, iscore, chitype):
        self.name = name
        self.p1 = p1        # first parent
        self.p2 = p2        # second parent
        self.p3 = p3        # third parent
        self.d1 = d1        # distance from first parent
        self.d2 = d2        # angle made with first and second parents
        self.d3 = d3        # torsion made withe first, second and third parents
        self.iscore = (iscore == '+' or iscore == 1)
        self.chitype = chitype

    def copy(self):
        return Atom(self.name, self.p1, self.p2, self.p3, self.d1,
                    self.d2, self.d3, self.iscore, self.chitype)

class Residue:
    """Class to hold residue information, adapted from LINUS"""
    def __init__(self, name, atoms=None):
        self.name = name
        if not atoms is None:
            self.atoms = atoms
        else:
            self.atoms = []
    def add_atom(self, atom):
        self.atoms.append(atom)

    def __getitem__(self, i):
        return self.atoms[i]

    def __len__(self):
        return len(self.atoms)

    def chiatom(self, chiname):
        for atom in self.atoms:
            if atom.chitype == chiname:
                return atom

    def copy(self):
        newRes = Residue(self.name)
        for atom in self:
            newRes.add_atom(atom.copy())
        return newRes

to_three_letter_code={'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
           'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
           'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
           'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
           'B': 'ASX', 'Z': 'GLX'}

def int2vec(v1, dis, v2, ang, v3, tors, sin=math.sin,\
            cos=math.cos, sqrt=math.sqrt):
    """int2vec(v1, dis, v2, ang, v3, tors)
        Compute coordinates for vector v, given its 'dis'tance
        from vector v1, 'ang'le v-v1-v2 and 'tors'ion v-v1-v2-v3.'
        Both ang and tors should be in degrees.
    """

    DEG2RAD = math.pi/180
    ang = ang*DEG2RAD
    sina = sin(ang)
    cosa = cos(ang)

    tors = (tors*DEG2RAD)
    sint = sina * sin(tors)
    cost = sina * cos(tors)

    x1, y1, z1 = v1
    x2, y2, z2 = v2
    x3, y3, z3 = v3

    u1x = x2 - x3
    u1y = y2 - y3
    u1z = z2 - z3
    d = 1.0 / sqrt(u1x*u1x + u1y*u1y + u1z*u1z)
    u1x = u1x*d
    u1y = u1y*d
    u1z = u1z*d

    u2x = x1 - x2
    u2y = y1 - y2
    u2z = z1 - z2;
    d = 1.0 / sqrt(u2x*u2x + u2y*u2y + u2z*u2z)
    u2x = u2x*d
    u2y = u2y*d
    u2z = u2z*d

    cosine = u1x*u2x + u1y*u2y + u1z*u2z

    if (abs(cosine) < 1.0):
        sine = 1.0/sqrt(1.0 - cosine*cosine)
    else:
        sine = 1.0/sqrt(cosine*cosine - 1.0)

    u3x = sine * (u1y*u2z - u1z*u2y)
    u3y = sine * (u1z*u2x - u1x*u2z)
    u3z = sine * (u1x*u2y - u1y*u2x)
    u4x = cost * (u3y*u2z - u3z*u2y)
    u4y = cost * (u3z*u2x - u3x*u2z)
    u4z = cost * (u3x*u2y - u3y*u2x)

    return [x1 + dis*(-u2x*cosa + u4x + u3x*sint),
            y1 + dis*(-u2y*cosa + u4y + u3y*sint),
            z1 + dis*(-u2z*cosa + u4z + u3z*sint)]

try:
    from vecmod import int2vec
except ImportError:
    pass


def seqtorib(name,sequence):
    """Converts a protein amino acid sequence into a LINUS-formated 'ribosome' file"""
    ribarray=[]
    header=['TITLE      %s'%(name),'DEFAULT PHI -120.0','DEFAULT PSI  120.0']
    ribarray.extend(header)
    for aa in sequence:
        if aa == 'P':
            ribarray.append('RES PRO PHI -70.0')
        else:
            ribarray.append('RES '+to_three_letter_code[aa])
    return "\n".join(ribarray)

def ribtoreslist(ribstring):
    allh=0
    residues=[]
    handle=StringIO.StringIO(ribstring)
    line=handle.readline()
    while line:
        if line.startswith('#'): continue
        line=line.strip().upper()
        if line.startswith('RES'):
            parts=line.split()
            name=parts[1]
            args=parts[2:]
            res = {'name':name}
            if args:
                for i in range(0,len(args),2):
                    res[args[i]]=float(args[i+1])
            residues.append(res)
        elif line.startswith('TITLE'):
            title=line
        elif line.startswith('DEF'):
            parts=line.split()
            ang=parts[1]
            val=float(parts[2])
            if ang == 'PHI': defphi = val
            if ang == 'PSI': defpsi = val
        elif line.startswith('ALLH'):
            allh=1
        line=handle.readline()
    handle.close()
 
    if defpsi:
        defpsi=defpsi+180.0
        if defpsi > 180: defpsi=defpsi-360.0

    for res in residues:
        if not res.has_key('PHI'):
            if defphi: res['PHI'] = defphi

        if not res.has_key('PSI'):
            if defpsi: res['PSI'] = defpsi
        else:
            val=res['PSI']+180.0
            if val>180:res['PSI']=val-360
            else: res['PSI']=val
    return title,allh,residues

def build_atom_from_rib_line(line):
    """Build an atom from a line in the parameter file"""
    args=line.upper().split()
    return Atom(args[0], int(args[4]), int(args[5]), int(args[6]),
                float(args[1]),float(args[2]), float(args[3]), args[7],
                args[8])

def read_parameters_from_file(paramfile):
    if type(paramfile)==type(''):
        paramdata=open(paramfile,'r')
    else: paramdata=paramfile
    line=paramdata.readline()
 
    residues={}

    while(line):
        line=line.strip()
        line=line.upper()
        if not line.startswith('NAME'):
            line=paramdata.readline()
            continue
        line=line.split()
        name=line[1]
        numat=int(line[-1])
        newres=Residue(name)
        for i in range(numat):
            line = paramdata.readline()
            newres.add_atom(build_atom_from_rib_line(line))
        residues[name]=newres
        line=paramdata.readline()
    return residues

def make_chain_from_db(residues,db):
    chain=[]
    for res in residues:
        name=res['name']
        if not db.has_key(name):
            raise TypeError,'Unknown Residue: %s'%name
        Res=db[name].copy()
        keys=res.keys()
        keys.remove('name')
        for key in keys:
            val=res[key]
            for atom in Res.atoms:
                if atom.chitype==key:
                    atom.d3=val
        chain.append(Res)
    return chain

def build_peptide(chain):
    def normalize_indices(chain):
        rp = map(len, chain)
        rp.insert(0, 0)
        rp = CUMSUM(rp)
        for i in range(len(chain)):
            res = chain[i]
            for atom in res:
                p = atom.p1
                if p < 0:
                    atom.p1 = rp[i-1] + abs(p) - 1
                else:
                    atom.p1 = rp[i] + p - 1
                p = atom.p2
                if p < 0:
                    atom.p2 = rp[i-1] + abs(p) - 1
                else:
                    atom.p2 = rp[i] + p - 1
                p = atom.p3
                if p < 0:
                    atom.p3 = rp[i-1] + abs(p) - 1
                else:
                    atom.p3 = rp[i] + p - 1
    normalize_indices(chain)
    dtor=math.pi/180
    res0=chain[0]
    res0[0].xyz=ZEROS(3,'d')
    a2=res0[1]
    res0[1].xyz=ARRAY((a2.d1,0,0),'d')
    a3 = res0[2]
    if a3.p1 == 1:
        x = a2.d1 - a3.d1*math.cos(dtor*a3.d2)
    else:
        x = a3.d1*math.cos(dtor*a3.d2)
    y = a3.d1 * math.sin(a3.d2*dtor)
    res0[2].xyz = ARRAY((x, y, 0.0), 'd')
    atomlist = []
    for res in chain:
        for atom in res.atoms:
            atomlist.append(atom)

    for i in xrange(3, len(atomlist)):
        a = atomlist[i]
        p1 = atomlist[a.p1].xyz
        p2 = atomlist[a.p2].xyz
        p3 = atomlist[a.p3].xyz
        a.xyz = int2vec(p1, a.d1, p2, a.d2, p3, a.d3)

def ribosome_to_coord(ribtup,paramfile):
    title,allh,residues=ribtup 
    db=read_parameters_from_file(paramfile)
    chain = make_chain_from_db(residues,db)
    build_peptide(chain)
    return chain_to_pdb(title,chain,allh)

def chain_to_pdb(title,chain,allh):
    f=StringIO.StringIO()
    if title.startswith('TITLE'):title=" ".join( title.split()[1:] )
    title=title.capitalize()

    fmt1 = 'ATOM  %5d  %-4s%-4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
    fmt2 = 'ATOM  %5d %-4s %-4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n'
    f.write('COMPND    %s\n' % title)
    seq = []
    nres = len(chain)
    line_num = 1
    for i in xrange(0, nres, 13):
        prefix = 'SEQRES%4d  %4d  ' % (line_num, nres)
        rec = [prefix]
        for j in range(i, min(i+13, nres)):
            rec.append(' '+chain[j].name)
        rec.append('\n')
        f.write(''.join(rec))
        line_num = line_num + 1
    i = 1
    j = 1
    for res in chain:
        res_name = res.name
        for atom in res:
            name = atom.name
            if name == 'DU': continue
            if atom.iscore or allh:
                if name[:1] in '123':
                    fmt = fmt2
                else:
                    fmt = fmt1
                f.write(fmt % ( i, name, res_name, j, atom.xyz[0], atom.xyz[1],
                                atom.xyz[2], 1.0, 0.0))
                i = i + 1
        j = j + 1
    f.write('END\n')
    f.seek(0)
    parser=PDB.PDBParser()
    protein=parser.get_structure(title,f)
    f.close()
    return protein
   
def seq_to_PDB(infile,format="Fasta",paramfile=None):
    if not paramfile:
        import urllib
        paramfile=urllib.urlopen('http://protlib.svn.sourceforge.net/viewvc/protlib/data/ribosome.dat')
        
    ribstrings={}
    structures={}

    extension=os.path.splitext(infile)[1]
    if extension.startswith(os.extsep): extension=extension[len(os.extsep):]
  
    if format: extension=format

    if type(infile)==type(''):
        try:
            infile=os.path.abspath(os.path.expanduser(infile))
            handle = open(infile, "rU")
        except IOError,err:
            print str(err),"require a valid input file"
    else: handle=infile
     
    record_dict = SeqIO.to_dict(SeqIO.parse(handle,extension))
    handle.close()
    for prot in record_dict.keys():
        ribstrings[prot]=seqtorib(prot,record_dict[prot].seq)
        structures[prot]=add_hn(ribosome_to_coord(ribtoreslist(ribstrings[prot]),paramfile))
    return structures

def add_hn(prot):
    def normalize(vec):
        vecsq=vec*vec
        return vec/numpy.sqrt(vecsq.sum())
    atomlist=[]
    for residue in prot.get_residues():
        for atom in residue:
            if atom.get_name()=="N": natom=atom
            if atom.get_name()=="CA": caatom=atom
            if atom.get_name()=="C": catom=atom
        atomlist.append([natom,caatom,catom])

    for resnum,residue in enumerate( prot.get_residues() ):
        if resnum==0 or residue.get_resname()=="PRO": continue
        catom=atomlist[resnum-1][2]
        natom=atomlist[resnum][0]
        caatom=atomlist[resnum][1]

        ncavec=caatom.get_coord()-natom.get_coord()
        ncvec=catom.get_coord()-natom.get_coord()
        ncavec=normalize(ncavec)
        ncvec=normalize(ncvec)
        sumvec=ncavec+ncvec
        sumvec=normalize(sumvec)
        NHvec=-sumvec*NHLength+natom.get_coord()
        nhatom=PDB.Atom.Atom('HN',NHvec,1,1," "," HN ",1)
        residue.add(nhatom)
    return prot

def usage():
    usagestring="""
 
    Usage: %s -h -i sequence_file -d data_file

              -h/--help        print this information
              -i/--input       input file with biopython SeqIO format as extension (e.g. 'fasta')
              -d/--data        residue information file 'ribosome.dat' from LINUS
    """ % (sys.argv[0])

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],"hi:d:",["help","input=","data="])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    infile="" 
    datafile=None

    for o,a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-i","--input"):
            infile=a
        elif o in ("-d","--data"):
            datafile=a
        else:
            assert False,"unhandled option"

    if not infile:
        print "Required file: sequence file"
        usage()
        sys.exit(2)
   
    result=seq_to_PDB(infile,datafile)
    for key in result.keys():
        io=PDB.PDBIO()
        io.set_structure(result[key])
        io.save(sys.stdout)

if __name__=="__main__":
    main()
