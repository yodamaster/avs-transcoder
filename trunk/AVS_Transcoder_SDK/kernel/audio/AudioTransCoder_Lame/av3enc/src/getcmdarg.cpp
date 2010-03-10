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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "getcmdarg.h"

static int isNum(char* opt)
{
    /* check if normal numeric expression: (+/-)xxxx.xxxx; (+/-).xxxx 
       or scientific numeric expression: (+/-)xxxx.xxxx(E/e)(+/-)xxxx
       or xxxx.xxxx% */    

    int     i;
    int     needSign;
    int     needPoint;
    int     needExp;
    int     needPer;
    
    i = 0;
    needSign = 1;
    needPoint = 1;
    needExp = 1;
    needPer = 0;
            
    while(opt[i] != '\0'){
        if((opt[i] == '+' || opt[i] == '-') && needSign)
            needSign = 0;
        else if(isdigit(opt[i])){
            needSign = 0;
            needPer = 1;
        }else if((opt[i] == 'E' || opt[i] == 'e') && needExp){
            needSign = 1;
            needExp = 0;
            needPer = 0;
        }else if(opt[i] == '.' && needPoint){
            needSign = 0;
            needPoint = 0;
            needPer = 0;    
        }else if(opt[i] == '%' && needPer){
            needSign = 0;
            needPoint = 0;
            needExp = 0;
            needPer = 0;            
        }else
            return 0;
        
        i++;
    }
    
    return 1;
}

static int isOpt(char* opt, char** optSet)
{
    /* check if 'opt' is an elegible option for switch indexed 'swIdx' 
       in case finite set option */
     
    if(opt == NULL)
        return 1;
    else
    {
        int i;
        
        i = 0;
        while(strcmp("\0", optSet[i])){
            if(!strcmp(optSet[i], opt))
                return 1;
            i++;
        }
        return 0;
    }
}

    
int parseCommandLine(cmd_option*         option,
                     const cmd_switch*   swts,
                     const int           argc, 
                     char*               argv[])
{
    int i;    
    int opIdx;
    
    int     numSwitch;
    char*   haveSet;
    
    numSwitch = 0;
    while(strcmp("\0", swts[numSwitch].swt)) numSwitch++;
    haveSet = (char *)malloc(numSwitch*sizeof(char));
    memset(haveSet, 0, numSwitch*sizeof(char));               
    
    i = 1;
    opIdx = 0;     
    while(i < argc){
        char mark;
        
        mark = argv[i][0];
        if(mark != '/' && mark != '-') {
            /* input files */
            option[opIdx].swIdx = numSwitch;
            option[opIdx].opt = argv[i++];
            opIdx++;
        }else{
            char* swt;
            char* opt;
            
            int j;
            int len;
            int s;
            int nst;
 
            /* go to the start of a new switch, all '/' or '-' between igored */
            swt = argv[i];  
            while(swt[0] == '/' || swt[0] == '-'){
                if(strlen(swt) == 1){
                    if(++i == argc){
                        /* come to the end of the cmd line */
                        printf("command line warning: \tno switch after '%c', '%c' ignore!\n", mark, mark);
                        return 0;
                    }else{
                        swt = argv[i];
                    }
                }else{
                    swt += 1;
                }
            }
            
            /* get switch idx and the start of the option for this switch 
               if it is contained in the current argment */
            
            s = -1;
            j = 0;
            nst = 1;
            opt = NULL;
            len = strlen(swt);            
            while(j < numSwitch){
                if(swts[j].optType == 0){
                    if(!strcmp(swts[j].swt, swt)){
                        s = j;
                        break;
                    }
                }else{
                    int slen;

                    slen = strlen(swts[j].swt);
                    if(!memcmp(swts[j].swt, swt, slen)){
                        s = j;
                        if(slen < len){
                            opt = swt + slen;
                            nst = 0;
                        }
                        break;
                    }
                }
                j++;
            }

            if(s < 0){
                /* no elegible switch found, ignore */
                printf("command line warning: \tunrecognized switch '%c%s', command line parsing restarts at the next argument!\n", mark, swt);
                i++;
                continue;
            }
            
            /* get the start of the option for this switch if it starts at a new argv */ 
            if(swts[s].optType){ 
                if(opt == NULL){
                    /* option starts at the next argv */
                    if(++i < argc)
                        opt = argv[i];
                    else if(swts[s].optType != 3)
                    {
                        /* come to end of cmd line while expecting an argument */
                        printf("command line warning: \t'%c%s' needs an argument, '%c%s' ignored!\n", mark, swts[s].swt, mark, swts[s].swt);
                        return opIdx;
                    }
                }
            }else if(opt != NULL){
                /* switch without option */
                printf("command line warning: \t'%c%s' does not need any argument, '%s' ignored!\n", mark,swts[s].swt, opt);
                opt = NULL;
            }
 
            /* check if argument valid */
            if(swts[s].optType == 1 && !isNum(opt)){
                /* no valid input found */
                printf("command line warning: \t'%c%s' need a numeric arguament, default value taken!\n", mark, swts[s].swt);
                if(!nst)
                    i++;
                continue;
            }
            if(swts[s].optType == 3 && !isOpt(opt, swts[s].optSet)){
                /* no valid input found */
                if(!nst)
                    printf("command line warning: \t'%s' is not an elegible argument for '%c%s', '%s' ignored!\n", opt, mark, swts[s].swt, opt);                 
                opt = NULL;
                i -= nst;
            }
            
            /* check if repeating setting the same nonaccumlative option */ 
            if(!swts[s].setType){
                if(haveSet[s])
                    printf("command line warning: \tcommand line parameter corresponding '%c%s' has already been set, \n\t\t\tthe new value will overwrite the old one!\n", mark, swts[s].swt);
                else
                {
                    int m;
                    
                    m = s;
                    while(!haveSet[m]){
                        haveSet[m] = 1;
                        m = swts[m].relative;
                    }
                }
            }  
            
            option[opIdx].swIdx = s;
            option[opIdx].opt = opt;
            opIdx++;
            i++;
        }
    }
    
    free(haveSet);
    return opIdx;
}

                                                                              