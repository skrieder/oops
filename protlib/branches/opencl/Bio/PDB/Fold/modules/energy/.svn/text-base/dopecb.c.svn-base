/**
 Protein Folding Library.
 * Copyright (C) 2009 Andres Colubri.
 * Contact: andres.colubri 'AT' gmail.com
 *
 * This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
 * Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
 * Homepage: http://ibd.uchicago.edu/
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * DOPE-CB statisical potential function.
 * Authors: Min-yi Shen (original creator of DOPE), Andres Colubri (protlib implementation), Glen Hocky (protlib2 port)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../core/memory.h"
#include "../../core/fileio.h"
#include "../../core/errorhandler.h"
#include "dopecb.h"

#define MAXBINNUM 29
#define MAXATCODE 7
#define MAXAACODE 19

IndexValue NUMBINS=MAXBINNUM+1;
FloatValue DISTANCECUTOFF=15.0;
FloatValue BINSIZE;

static const char * aacodetoonelet[21]={" ","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};
static const char * dopecodetolet[8]={"N","HN","CA","CB","HA2","HA1","C","O"};
static const char * atcodetolet[10]={" ","N","CA","C","O","HN","HA1","CAP","CB","HA2"};
static const IndexValue convtodope[10]={-1,0,2,6,7,1,5,-1,3,4};

struct AtomData *atoms;          // Array of atom data.
FloatValue *coords;              // Array of atom coordinates.
FloatValue *distances;           // Array of atom pair-wise distances.
FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
BoolValue *atommask;             // Array of atom masking values.
BoolValue *distmask;             // Array of distance masking values.  
IndexValue natoms;               // Total number of atoms.

struct ResidueData *residues;    // Array of residue data.
IndexValue nres;                 // Total number of residues.

struct ChainData *chains;        // Array of chain data.
IndexValue nchains;              // Total number of chains.

IndexValue natoms2,dopearraysize;

FloatValue * dopeenergies;
FloatValue * pairenergies;

IndexValue * pairrows;
IndexValue * activepairs;
IndexValue nactivepairs;
IndexValue * updatepairs;
IndexValue nupdatepairs;

BoolValue debug;

ErrorCode print_pair_data(IndexValue at1,IndexValue at2);
ErrorCode print_pair_data(IndexValue at1,IndexValue at2){
    IndexValue res1=atoms[at1].res;
    IndexValue res2=atoms[at2].res;
    IndexValue at1type=atoms[at1].type;
    IndexValue at2type=atoms[at2].type;
    IndexValue res1type=residues[res1].type;
    IndexValue res2type=residues[res2].type;
    IndexValue dopecode1=atoms[at1].niprop;
    IndexValue dopecode2=atoms[at2].niprop;
    //for full detail, uncomment
    //printf("%s(%i) %s(%i) %s(%i) %s(%i) %i %i",aacodetoonelet[res1type],res1type,dopecodetolet[dopecode1],dopecode1,aacodetoonelet[res2type],res2type,dopecodetolet[dopecode2],dopecode2,res1+1,res2+1);
    printf("%s %s %s %s %i %i",aacodetoonelet[res1type],dopecodetolet[dopecode1],aacodetoonelet[res2type],dopecodetolet[dopecode2],res1+1,res2+1);

    return NO_ERROR;
}

ErrorCode set_system_for_dopecb(struct AtomData *_atoms, FloatValue *_coords,
                                                         FloatValue *_distances, FloatValue *_invdistances,
                                                         BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms,
                                struct ResidueData *_residues, IndexValue _nres,
                                struct ChainData *_chains, IndexValue _nchains)
{
    atoms = _atoms;
    coords = _coords;
    distances = _distances;
    invdistances = _invdistances;
    atommask = _atommask;
    distmask = _distmask;
    natoms = _natoms;
    residues = _residues;
    nres = _nres;
    chains = _chains;
    nchains = _nchains;

    return NO_ERROR;
}

ErrorCode load_dopecb_parfile(const char* parfilename){
  /*IO ideas from http://www.mrx.net/c/source.html*/

  FILE *inputfilehandle=fopen(parfilename,"r");

  if (inputfilehandle == NULL) {
      print_error(ITEM_NOT_FOUND_ERROR,"load_dopecb_parfile", 0);
  }

  dopearraysize=convert_to_dopecb_pos(MAXAACODE,MAXATCODE,MAXAACODE,MAXATCODE,MAXBINNUM);

  create_float_array(&dopeenergies,dopearraysize);

#define CR 13            /* Decimal code of Carriage Return char */
#define LF 10            /* Decimal code of Line Feed char */
#define EOF_MARKER 26    /* Decimal code of DOS end-of-file marker */
#define MAX_REC_LEN 1024 /* Maximum size of input buffer */

  /*start of T_fread*/
  IntValue  isNewline;              /* Boolean indicating we've read a CR or LF */
  IntValue  lFileLen;               /* Length of file */
  IntValue  lIndex;                 /* Index into cThisLine array */
  IntValue  lLineCount;             /* Current line number */
  IntValue  lLineLen;               /* Current line length */
  IntValue  lStartPos;              /* Offset of start of current line */
  IntValue  lTotalChars;            /* Total characters read */
  char cThisLine[MAX_REC_LEN]; /* Contents of current line */
  char *cFile;                  /* Dynamically allocated buffer (entire file) */
  char *cThisPtr;               /* Pointer to current position in cFile */

  fseek(inputfilehandle, 0L, SEEK_END);  /* Position to end of file */
  lFileLen = ftell(inputfilehandle);     /* Get file length */
  rewind(inputfilehandle);               /* Back to start of file */

  cFile = calloc(lFileLen + 1, sizeof(char));

  if(cFile == NULL )
  {
    printf("\nInsufficient memory to read file.\n");
    return 0;
  }

  fread(cFile, lFileLen, 1, inputfilehandle); /* Read the entire file into cFile */

  lLineCount  = 0L;
  lTotalChars = 0L;

  cThisPtr    = cFile;              /* Point to beginning of array */

  while (*cThisPtr)                 /* Read until reaching null char */
  {
    lIndex    = 0L;                 /* Reset counters and flags */
    isNewline = 0;
    lStartPos = lTotalChars;

    while (*cThisPtr)               /* Read until reaching null char */
    {
      if (!isNewline)               /* Haven't read a CR or LF yet */
      {
        if (*cThisPtr == CR || *cThisPtr == LF) /* This char IS a CR or LF */
          isNewline = 1;                        /* Set flag */
      }

      else if (*cThisPtr != CR && *cThisPtr != LF) /* Already found CR or LF */
        break;                                     /* Done with line */

      cThisLine[lIndex++] = *cThisPtr++; /* Add char to output and increment */
      ++lTotalChars;

    } /* end while (*cThisPtr) */

    cThisLine[lIndex] = '\0';     /* Terminate the string */
    ++lLineCount;                 /* Increment the line counter */
    lLineLen = strlen(cThisLine); /* Get length of line */

    /* Print the detail for this line */
    //printf("%s",cThisLine);
    //PrintLine(cThisLine, lLineCount, lLineLen, lStartPos, NULL, 0);

    char * pch;
    IndexValue count=0;
    IntValue indices[4];
    IntValue currentindex;
    pch = strtok (cThisLine," ");
    while (pch != NULL)
    {
      if (count<4) {
          indices[count]=atoi(pch);
          //printf("%i ",indices[count]);
          }
      else {
          FloatValue currentenergy=atof(pch);
          //printf(" %f ",currentenergy);
          currentindex=convert_to_dopecb_row(indices[0],indices[1],indices[2],indices[3])+ (count-4);
          dopeenergies[currentindex]=currentenergy;
      }
      pch=strtok(NULL," ");
      count++;
    }
    //printf(" %i\n",count);
  } /* end while (cThisPtr <= cEndPtr) */

  //printf("Length of file array=%#x (dec %d)\n", strlen(cFile), strlen(cFile));

  free(cFile);
/*end of T_fread*/

  fclose(inputfilehandle);

  return NO_ERROR;
}

ErrorCode init_dopecb(const char* cfg)
{
    debug=0;
    nactivepairs=0;
    natoms2=natoms*natoms;

    create_float_array(&pairenergies,natoms2);

    create_index_array(&pairrows,natoms2);
    create_index_array(&updatepairs,natoms2);
    create_index_array(&activepairs,natoms2);

    BINSIZE=DISTANCECUTOFF/NUMBINS;

    printf("        Initializing DOPE function...\n");
    printf("              Loading Parameter file...\n");

    char *linehead = NULL;
    char *parfilename = NULL;
    if (open_input_text_file(cfg) == NO_ERROR)
    {
        while (read_input_line())
        {
            get_input_strval(&linehead, "=", 0);

            char firstchar=linehead[0];
            if (firstchar=='#') {
                free(linehead);
                linehead = NULL;
                continue;
            }
            if (strcmp(linehead, "PARFILE") == 0)
            {
                get_input_strval(&parfilename, "=,", 1);
            }
            else if (strcmp(linehead, "DEBUG") == 0)
            {
                get_input_ival(&debug, "=,", 1);
            }
            else
            {
                print_error(UNKNOWN_PARAMETER_ERROR, "init_debug", 0);
                printf("...line: %s\n",linehead);
            }
            free(linehead);
            linehead = NULL; // Make sure of set this string pointer to NULL, beacause it will be used again when reading the next line. 
        }
        close_input_file();
    }
    load_dopecb_parfile(parfilename);
    free(parfilename);

    //check to make sure dope is loaded correctly
    IndexValue i,j,k,l,m;
    if (debug){
        for (i=0;i<MAXAACODE;i++)
            for(j=0;j<MAXATCODE;j++)
                for(k=0;k<MAXAACODE;k++)
                    for(l=0;l<MAXATCODE;l++){
                        //printf("%i %i %i %i ",i,j,k,l);
                        printf("%s(%i) %s(%i) %s(%i) %s(%i)",aacodetoonelet[i+1],i,dopecodetolet[j],j,aacodetoonelet[k+1],k,dopecodetolet[l],l);
                        for(m=0;m<=MAXBINNUM;m++){
                            printf(" %4.2f ",dopeenergies[convert_to_dopecb_row(i,j,k,l)+m]);
                        }
                        printf("\n");
                    }
    }
    IndexValue resi,typei,restypei,resj,typej,restypej;
    for (i=0;i<natoms;i++){
        IndexValue dopeattypei=convtodope[atoms[i].type];
        typei=atoms[i].niprop=dopeattypei;
        resi=atoms[i].res; // starts from 0
        restypei=residues[resi].type-1; //starts from 1
        //atoms[i].niprop=0;
        //if (debug) printf("%i) Atom: %i, Type: %i,Residue type: %i\n",resi,i,typei,restypei);
        for (j=0;j<natoms;j++){
            IndexValue dopeattypej=convtodope[atoms[j].type];
            typej=atoms[j].niprop=dopeattypej;
            resj=atoms[j].res; // starts from 0
            restypej=residues[resj].type-1; //starts from 1
            pairrows[natoms*i+j]=pairrows[natoms*j+i]=convert_to_dopecb_row(restypei,typei,restypej,typej);
            //printf("Dope Row: %i,idx1: %i, idx2: %i\n",row,idx1,idx2);
            if (resj>resi+0) { //change this to reflect ii+n only
                activepairs[2*nactivepairs]=i;
                activepairs[2*nactivepairs+1]=j;
                nactivepairs+=1;
            }
            pairenergies[natoms*i+j]=0;
        }
    }
    printf("        Done.\n");
    return NO_ERROR;
}

ErrorCode finish_dopecb(void)
{
    delete_float_array(&pairenergies,natoms2);
    delete_index_array(&pairrows,natoms2);
    delete_index_array(&updatepairs,natoms2);
    delete_index_array(&activepairs,natoms2);

    delete_float_array(&dopeenergies,dopearraysize);
    return NO_ERROR;
}

ErrorCode pre_dopecb_calc(void)
{
    IndexValue i,j;
    nupdatepairs=0;
    for (i=0;i<natoms;i++)
        for (j=i;j<natoms;j++)
        {
            if (distmask[distance_index(i,j,natoms)]){
                updatepairs[2*nupdatepairs]=i;
                updatepairs[2*nupdatepairs+1]=j;
                nupdatepairs+=1;
            }
        }
    if (debug) printf("Update pairs: %i\n",nupdatepairs);
    return NO_ERROR;
}

ErrorCode post_dopecb_calc(void)
{
    return NO_ERROR;
}

ErrorCode calc_dopecb_energy(FloatValue *energy)
{
    FloatValue totalenergy=0.0;
    FloatValue r1,f1,f2;
    IndexValue i,pairidx,pairidxt,ri,idx1,idx2;

    for (i=0;i<nupdatepairs;i++){
        idx1=updatepairs[2*i];
        idx2=updatepairs[2*i+1];

        pairidx=natoms*idx1+idx2;
        pairidxt=natoms*idx2+idx1;

        IndexValue doperow=pairrows[pairidx];

        FloatValue dist=distances[distance_index(idx1,idx2,natoms)];
        if (0< dist && dist<DISTANCECUTOFF){
            r1 = dist/BINSIZE - 0.5;
            ri = floor(r1);

            f1 = r1 - ri;
            f2 = 1.0 - f1;

            f1=(ri==MAXBINNUM)?0:f1;           

            if (ri>0) {
                pairenergies[pairidx]=f2*dopeenergies[doperow+ri] + f1*dopeenergies[doperow+ri+1];
            }
            else pairenergies[pairidx] = dopeenergies[doperow];
        }
        else pairenergies[pairidx]=0.0;
/*
        printf("%i %i ",idx1,idx2);
        print_pair_data(idx1,idx2);
        printf(" %f %i %i %f %f %f\n",dist,ri,ri+1,f1,f2,pairenergies[pairidx]);//,dopeenergies[doperow+ri],pairenergies[pairidx]);

*/
        pairenergies[pairidxt]=pairenergies[pairidx];
    }
    //for (i=0;i<natoms*natoms;i++) printf("%i %f\n",i,pairenergies[i]);

    for (i=0;i<nactivepairs;i++){
        idx1=activepairs[2*i];
        idx2=activepairs[2*i+1];
        pairidx=natoms*idx1+idx2;

//        printf("%i %i ",idx1,idx2);
/*
        print_pair_data(idx1,idx2);
        printf(" %5.3f %4.2f\n",distances[distance_index(idx1,idx2,natoms)],pairenergies[pairidx]);
*/

        totalenergy+=pairenergies[pairidx];
    }
    *energy=totalenergy;
    return NO_ERROR;
}

ErrorCode calc_dopecb_gradient(FloatValue *gradient)
{
    return NO_ERROR;
}
