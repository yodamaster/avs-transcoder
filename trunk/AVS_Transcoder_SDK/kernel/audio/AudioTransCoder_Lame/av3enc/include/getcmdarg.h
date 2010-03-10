/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
* edited by
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the users without any
* license fee or royalty on an "as is" basis. The AVS disclaims
* any and all warranties, whether express, implied, or statutory,
* including any implied warranties of merchantability or of fitness
* for a particular purpose. In no event shall the contributors or 
* the AVS be liable for any incidental, punitive, or consequential
* damages of any kind whatsoever arising from the use of this program.
*
* This disclaimer of warranty extends to the user of this program
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The AVS does not represent or warrant that the program furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy is available from the AVS Web site at
* http://www.avs.org.cn
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

#ifndef __GETCMDARG_H
#define __GETCMDARG_H

#ifdef __cplusplus
extern "C"
{
#endif 

typedef struct {
    char*           swt;        /* switch word*/   
    int             optType;    /* 0, without option;  1, with a numeric option 
                                 2, with a string option; 3, finite set option */
    char**          optSet;     /* in case finite set option */
    int             setType;    /* 0, new set overwrites old one
                                 1, all sets are accummulated */
    int             relative;   /* the idx of the next switch affects the same cmd_params members */                                    
} cmd_switch;

typedef struct {        
    int     swIdx;      /* the idx of a switch */
    char*   opt;        /* the pointer of an option for the corresponding switch */
} cmd_option;

/* command line parsing strategy:
 * 0. a string without a '/' or '-' immediately preceding it or the first char is 
      not '/' or '-' is taken as ordinary string whatever it looks like.
 * 1. a switch has at most one argument.
 * 2. a string immediately follows '/' or '-' is taken as a switch; a string excluding
      the first char is also a switch if the one is '/' or '-'.
 * 3. the part of a string (if not vacant) immediately follows a switch is taken as
      as the argument of the switch if the switch expects an argument.
 * 4. a switch expecting an argument but followed by an invalid argument is ignored
      with warning. 
 * 5. an argument following a switch not expecting any argument is ignored with warning.
 * 6. an unrecognized switch is ignored with warning
 * 7. only the last seting of a cmd_params member except inFile and outFile is active. 
      the preceding setting (if any) is overwritten by the current one with warning.
 * 8. command line input about Infile and outFile are accumlated.
 * 9. mutiple '/' or '-' in a row simply ignored.
 
 * NOTE: TO AVOID SWITCH WORD AMBIGUITY, A SWITCH EXPECTING AN ARGUMENTS SHOULD NOT HAVE 
         ITS SWITCH WORD CONCINCIDE WITH THE LEADING CHARACTERS OF ANY OTHER SWITCHES.
         (i.e. "abc" and "abcd" should not appear in the same cmd_switch set if "abc" is
          attached to a switch expecting an argument)
 */

/* parse 'argc' number of argv to set 'param' according to 'swts'.
 * parsed options are stored in 'option'. return the number(>=0) of all active cmd line inputs
 * NOTE: the caller roution has the responsibility to ensure there is enough space
 *       in 'option' for storage.
 */ 
int parseCommandLine(cmd_option*         option,   /*out*/
                     const cmd_switch*   swts,     /*in*/
                     const int           argc,     /*in*/
                     char*               argv[]);  /*in*/
                            
#ifdef __cplusplus
}
#endif 

#endif 
