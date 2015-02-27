/****************************************************************\
*                                                                *
*  fastaoverlap: generate overlapping sequences                  *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2009.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include "globals.h"
#include "argument.h"
#include "fastadb.h"

FILE *file;

typedef struct {
    gint chunk_size;
    gint jump_size;
} FastaOverlap_Info;

static gboolean fasta_overlap_traverse_func(FastaDB_Seq *fdbs,
                                            gpointer user_data){
    FastaOverlap_Info *foi = user_data;
    gint region_start = 0, region_end = foi->chunk_size;
    Sequence *nfe;
    if(region_end > fdbs->seq->len)
        region_end = fdbs->seq->len;
    nfe = Sequence_subseq(fdbs->seq, region_start,
                          region_end - region_start);
    Sequence_print_fasta(nfe, file, FALSE);
    Sequence_destroy(nfe);
    while(region_end < fdbs->seq->len){
        region_start += foi->jump_size;
        if(region_start > fdbs->seq->len)
            break;
        region_end += foi->jump_size;
        if(region_end > fdbs->seq->len)
            region_end = fdbs->seq->len;
        nfe = Sequence_subseq(fdbs->seq, region_start,
                              region_end - region_start);
        Sequence_print_fasta(nfe, file, FALSE);
        Sequence_destroy(nfe);
        }
    return FALSE;
    }

int Argument_main(Argument *arg){
    FastaDB *fdb;
    ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *outputFile;
    FastaOverlap_Info foi;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'c', "chunk", NULL,
        "Maximum size of each chunk", "100000",
        Argument_parse_int, &foi.chunk_size);
    ArgumentSet_add_option(as, 'j', "jump", NULL,
        "Jump between each chunk", "90000",
        Argument_parse_int, &foi.jump_size);
    ArgumentSet_add_option(as, 'O', "output", "path",
        "Specify the output file", "stdout",
        Argument_parse_string, &outputFile);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaoverlap",
        "Generate overlapping fasta sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);

    if (g_strcmp0(outputFile, "stdout") != 0) {
        fprintf(stdout, "Writing output to %s\n", outputFile);
        file = fopen(outputFile, "w");
    } else {
        file = stdout;
    }
    if (file == NULL) {
        fprintf(stderr, "Could not create output file '%s'\n", outputFile);
        exit(-1);
    }

    fdb = FastaDB_open(query_path, NULL);
    if(foi.chunk_size < 1)
        g_error("Chunk size (%d) must be greater than zero",
                foi.chunk_size);
    if(foi.jump_size < 1)
        g_error("Jump size (%d) must be greater than zero",
                foi.jump_size);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_overlap_traverse_func, &foi);
    FastaDB_close(fdb);
    return 0;
    }
