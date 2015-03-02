/****************************************************************\
*                                                                *
*  Library for word-neighbourhood generation                     *
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

#include <string.h>
#include "wordhood.h"
#include "argument.h"

static gboolean wordhood_test_func(gchar *word, gint score,
                                   gpointer user_data){
    gint *count = user_data;
    g_message("id [%d] word [%s] score [%d]", ++(*count), word, score);
    return FALSE;
    }

int Argument_main(Argument *arg){
    gchar *seq = "AAACCCGGGTTT";
    Submat *s = Submat_create("nucleic");
    CodonSubmat *cs = CodonSubmat_create();
    WordHood_Alphabet *wha;
    WordHood *wh;
    gint count;
    /**/
    g_message("using nucleic submat");
    count = 0;
    wha = WordHood_Alphabet_create_from_Submat("ACGT", "ACGT",
                                               s, FALSE);
    wh = WordHood_create(wha, 9, TRUE);
    WordHood_info(wh);
    WordHood_traverse(wh, wordhood_test_func, seq, strlen(seq), &count);
    WordHood_destroy(wh);
    WordHood_Alphabet_destroy(wha);
    /**/
    g_message("using codon submat");
    count = 0;
    wha = WordHood_Alphabet_create_from_CodonSubmat(cs, FALSE);
    wh = WordHood_create(wha, 9, TRUE);
    WordHood_info(wh);
    g_message("begin trav");
    WordHood_traverse(wh, wordhood_test_func, seq, strlen(seq), &count);
    g_message("done trav");
    WordHood_destroy(wh);
    WordHood_Alphabet_destroy(wha);
    /**/
    Submat_destroy(s);
    CodonSubmat_destroy(cs);
    return 0;
    }

