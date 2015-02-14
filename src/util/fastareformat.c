/****************************************************************\
*                                                                *
*  fastalength : a utility to reformat fasta format sequences    *
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

static gboolean fasta_reformat_traverse_func(FastaDB_Seq *fdbs,
                                           gpointer user_data){
    FastaDB_Seq_print(fdbs, file, FastaDB_Mask_ID
                                   |FastaDB_Mask_DEF
                                   |FastaDB_Mask_SEQ);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *outputFile;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as_input, 'O', "output", "path",
        "Specify the output file", "stdout",
        Argument_parse_string, &outputFile);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastareformat",
        "A utility to reformat fasta sequence files\n"
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
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_reformat_traverse_func, NULL);
    FastaDB_close(fdb);
    return 0;
    }
