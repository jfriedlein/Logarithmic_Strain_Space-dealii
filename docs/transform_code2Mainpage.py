#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 18:50:30 2019

@author: jfriedlein
"""
# useful commands to clear the workspace and clear the console (MATLAB: clear,clc) from https://stackoverflow.com/questions/54943710/code-to-clear-console-and-variables-in-spyder
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#
# Transform normal code into Doxygen mainpage
#
    
# set the paths to ...
dirMainpage = "./mainpage-framework.h" # ... the framework file containing e.g. the introduction text
#dirCode = "../hybrid_solver.h"  # ... the code that will be placed into the mainpage between the keywords "\code" and "\endcode"
mergedFile = "../mainpage.h"  # ... the output file containing the text from the framework and the code from the dirCode-file

do_replace = False
with open(mergedFile, "w") as outputFile: # write into the output file line by line
    with open(dirMainpage, "r") as mainpageDoc: # read the framework file line by line until ...
        for line in mainpageDoc:
            if ( "\code" in line ):  # ... you find this keywork, then place the code after this keyword until ...
                # First we extract the name of the code file that should be placed between \code and \endcode
                next_line = next(mainpageDoc)
                dirCode = "../" + next_line[0:len(next_line)-1] # The last part just removes the \n at the end of the line
                # modify the code and write it into the mainpage file:
                code = True # True: the current line contains code;
                            # False: the current line contains a comment that shall be outputed as normal text
                firstLine = True
                with open(dirCode, "r") as codeDoc:  # read the code line by line
                    for lineCode in codeDoc:
                        # read the first few characters from this line that aren't blank
                        first_nonspace_index = len(lineCode) - len(lineCode.lstrip())
                        first_N_nonspace_letters = lineCode[first_nonspace_index:first_nonspace_index+3]
                        codeLast = code # stores whether the last line was code or comment
                        
                        # instead of /// we indicate a completely commented out line by use no space after //,
                        #  e.g. "//std::cout ..." is commented out but not a comment like "// This is a comment"
                        #if ( first_N_nonspace_letters == "///" ): # this indicates some completely outcommented code, so we keep the comment and don't create a text from this
                        #    lineCode = lineCode[1:len(lineCode)] # on
                        #    code = True
                        # ToDo-optimize: The following means that we only identify a line as comment when it exactly starts with "// ",
                        # with the emphasis on the blank space. When you use a tab direclty after the "//" it won't be identified as a comment.
                        # So you need to always use "// ...", where the dots can be as many tabs as you like, after you typed the blank space.
                        if ( first_N_nonspace_letters == "// " ):
                            if ( len(lineCode)==3 ): # this would be an empty line marked as comment
                                lineCode = "\n"
                            else:
                                lineCode = lineCode[(first_nonspace_index+3):len(lineCode)]
                            code = False
                        elif ( first_N_nonspace_letters == "//\n" ): # when it is an entirely empty comment line (without blank after "//")
                            lineCode = "\n"
                            code = False
                        elif ( first_N_nonspace_letters == "" ): # line is an empty line
                            lineCode =" \n"
                            code = codeLast     # save the type of the last line and do this until you find code or comment
                        else: # the line is actually code
                            # lineCode is used unchanged
                            code = True
                         
                        # At first, we have to either create a code-block for a line containing code or we start with normal text for comment
                        if ( firstLine == True ):
                            if ( code == True ): # if the first line is a code then we have to print the "\code"-command
                                codeLast = False # to trigger the corresponding "elif" below
                            elif ( code == False ):
                                codeLast = False # so we will just print the line as it is
                            firstLine = False
                                                
                        if ( code==codeLast ): # the current line is of the same type as the last, so we can remain in the currently active block (either code embraced in \code...\endcode or a text)
                            outputFile.write(lineCode)
                        elif ( codeLast==False and code==True ): # the text paragraph ended and now we found a line of code ...
                            outputFile.write("\code\n")  # ... so we start a code block
                            outputFile.write(lineCode)
                        elif ( codeLast == True and code == False ): # our last line was code and now we found a text ...
                            outputFile.write("\endcode\n") # ... so we close the current code block and start a normal text
                            outputFile.write(lineCode)
                            
                    if ( code == True ): # ensure that the blocks are closed at the end
                        outputFile.write("\endcode\n")
                do_replace = True
        
            # if we don't add the text from the code into the outputFile,
            #  we want to continue copying everything from the framework file into the output
            if ( do_replace == False ): 
                outputFile.write(line)
        
            # ... you find the "\endcode"-command in the source mainpage then stop replacing things,
            #     but only for the next line because we don't want to add this \endcode-command to the outputFile
            if ( do_replace == True and ("\endcode" in line) ): 
                do_replace = False
                

            
