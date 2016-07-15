#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys



def write_t_file(out_file, a):
    with open(out_file, 'w') as h:
        if type(a) is list:
            for i in range(len(a)):
                if type(a[i]) is list:
                    if len(a[i])>1:
                        for j in range(len(a[i])-1):
                            h.write(str(a[i][j])+"\t")
                        j = j+1
                        h.write(str(a[i][j])+"\n")
                    else:
                        if len(a[i]) == 1:
                            h.write(str(a[i][0])+"\n")
                        else:
                            h.write("\n")
                else:
                    h.write(str(a[i])+"\n")
        else:
            h.write(str(a)+"\n")


#write_file("result.txt", [["1"], [1,2]])
