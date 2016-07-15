#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

def parse_dist(in_file):
    result = []
    distribution = {}
    name = []
    with open(in_file, 'r') as f:
        flag = 0
        for aline in f.readlines():
            line = aline.strip()
            dist = line.split('\t')
            if len(dist[0].strip()) > 0:
                if len(dist) == 1:
                    if flag == 0:
                        name.append(line)
                        flag = 1
                        t_name = line
                    else:
                        distribution[t_name] = 'null'
                        name.append(line)
                        flag = 1
                        t_name = line
                else:
                    distri = []
                    for i in range(0, len(dist)):
                        distri.append(dist[i].strip())
                    distribution[t_name] = distri
                    flag = 0
    result.append(name)
    result.append(distribution)
    return result
                
                







        





